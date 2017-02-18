#!/usr/bin/env python2.7
from __future__ import print_function

import argparse
import os
import multiprocessing
import subprocess
import sys
import textwrap
import tarfile
from urlparse import urlparse
from contextlib import closing

import yaml
from bd2k.util.files import mkdir_p
from bd2k.util.processes import which
from toil.job import Job
from toil.lib.docker import dockerCall
from toil_lib import require, UserError
from toil_lib.files import tarball_files, copy_files
from toil_lib.jobs import map_job
from toil_lib.urls import download_url, s3am_upload

# formats
SCHEMES = ('http', 'file', 's3', 'ftp')
FILE_TYPES = ('fq', 'tar', 'bam')
PAIRED_TYPES = ('single', 'paired')
FQ_FORMATS = ('.fq', '.fastq', '.fq.gz', '.fastq.gz')
TAR_FORMATS = ('.tar', '.tar.gz')
BAM_FORMATS = ('.bam')
ZIP_FORMATS = ('.fq.gz', '.fastq.gz', '.tar', '.tar.gz')
PAIRED_FRONT_SUFFIXES = ("_1", "R1")
PAIRED_BACK_SUFFIXES = ("_2", "R2")

# filenames
DEFAULT_CONFIG_NAME = 'config-toil-topmed.yaml'
DEFAULT_MANIFEST_NAME = 'manifest-toil-topmed.tsv'

# docker images
#todo unlatest this
DOCKER_SAMTOOLS = "quay.io/ucsc_cgl/samtools:latest"

# todo temporary for development
trevdev = True

# Pipeline specific functions
def parse_samples(path_to_manifest):
    """
    Parses samples, specified in either a manifest or listed with --samples

    :param str path_to_manifest: Path to configuration file
    :return: Samples and their attributes as defined in the manifest
    :rtype: list[list]
    """

    samples = []
    with open(path_to_manifest, 'r') as f:
        for line in f.readlines():
            if line.isspace() or line.startswith('#'):
                continue
            sample = line.strip().split('\t')

            # validate structure
            if len(sample) != 4:
                raise UserError('Bad manifest format! Expected 4 tab separated columns, got: {}'.format(sample))

            # extract sample parts
            uuid, type, paired, given_url = sample
            url = None

            # validate contents
            if type not in FILE_TYPES:
                raise UserError("Incorrect file type '{}'.  Expected: {}".format(type, FILE_TYPES))
            if paired not in PAIRED_TYPES:
                raise UserError("Incorrect paired type '{}'.  Expected: {}".format(type, PAIRED_TYPES))

            # if URL is a local folder
            if urlparse(given_url).scheme == '':
                if type == 'bam':
                    valid_types = BAM_FORMATS
                elif type == 'fq':
                    valid_types = FQ_FORMATS
                else:
                    raise UserError("Files in local folder must be of type: ['bam', 'fq']")
                url = ['file://' + os.path.join(given_url, x) for x in filter(lambda x: x.endswith(valid_types), os.listdir(given_url))]
                require(len(url) > 0, "Found no files in local folder: '{}'".format(given_url))

            # If url is a tarball
            elif type == 'tar':
                url = given_url.split(',')
                require(len(url) == 1,
                        'tar URL "{}" not valid. "tar" samples require a single url.'.format(given_url))
                require(given_url.endswith(TAR_FORMATS),
                        'tar URL "{}" not valid. Approved extensions: {}'.format(given_url, TAR_FORMATS))
                require(urlparse(given_url).scheme in SCHEMES,
                        'tar URL "{}" not valid. Schemes: {}'.format(given_url, SCHEMES))

            # If URL is a fastq or series of fastqs
            elif type == 'fq':
                url = given_url.split(',')
                [require(u.endswith(FQ_FORMATS),
                         'fq URL "{}" not valid. Approved extensions: {}'.format(u, FQ_FORMATS)) for u in url]
                [require(urlparse(u).scheme in SCHEMES,
                         'fq URL "{}" not valid. Schemes: {}'.format(u, SCHEMES)) for u in url]

            # If URL is a bam or series of bams
            elif type == 'bam':
                url = given_url.split(',')
                [require(u.endswith(BAM_FORMATS),
                         'bam URL "{}" not valid. Approved extensions: {}'.format(u, BAM_FORMATS)) for u in url]
                [require(urlparse(u).scheme in SCHEMES,
                         'bam URL "{}" not valid. Schemes: {}'.format(u, SCHEMES)) for u in url]

            # sanity check
            else:
                raise UserError('PROGRAMMER ERROR: unexpected type "{}"'.format(type))

            # validate paired samples
            if paired == 'paired' and type == 'fq':
                r1 = len(_files_with_suffix(url, "R1"))
                r2 = len(_files_with_suffix(url, "R2"))
                u1 = len(_files_with_suffix(url, "_1"))
                u2 = len(_files_with_suffix(url, "_2"))
                require(r1 == r2, 'URL "{}" not valid for "paired" type: inconsistent files ending with "R1" and "R2"')
                require(u1 == u2, 'URL "{}" not valid for "paired" type: inconsistent files ending with "_1" and "_2"')
                require(r1 + r2 + u1 + u2 > 0,
                        'URL "{}" not valid for "paired" type: no files ending with "R1", "R2", "_1", or "_2"')

            sample = [uuid, type, paired, url]
            samples.append(sample)
    return samples


def _files_with_suffix(urls, suffix):
    """
    Returns filtered list of files in 'urls' which end with the suffix (disregarding file type)
    :param urls: list of fully-qualified file names (with or without scheme)
    :param suffix: string (or string tuple) to filter files on
    :return: list of strings
    """
    return filter(lambda u: os.path.splitext(os.path.basename(urlparse(u).path))[0].endswith(suffix), urls)

def prepare_input(job, sample, config):
    """
    input: list of files in .bam or .fq or tarballs
    for each input file:
        if tarball:
            unzip
        if bam
            picard sam2fq
        ensure paird or not paired (as specified in manifest)
    cat all fq's together
        (this will be a big file, is this an issue?)
    output: one (paired potentially) big fastq file
    """

    # prepare
    config = argparse.Namespace(**vars(config))
    config.cores = min(config.maxCores, multiprocessing.cpu_count())
    uuid, file_type, paired, urls = sample
    is_paired = paired == "paired"
    config.uuid = uuid
    work_dir = job.fileStore.getLocalTempDir()
    fastq_location = os.path.join(work_dir, "fastq")
    os.mkdir(fastq_location)
    bam_location = os.path.join(work_dir, "bam")
    os.mkdir(bam_location)

    # first, get all files
    for url in urls:
        #need to untar them (for .tar, .tar.gz, .fq.gz, .fastq.gz)
        if url.endswith(TAR_FORMATS):
            #todo: verify this works with all the tar formats (including fq stuff)
            tar_path = os.path.join(work_dir, os.path.basename(url))
            download_url(job, url=url, work_dir=work_dir)
            subprocess.check_call(['tar', '-xvf', tar_path, '-C', fastq_location])
            os.remove(tar_path)
        #need to extract to fq
        elif file_type == 'bam':
            # download bam
            download_url(job, url=url, work_dir=bam_location)
            # prep for extraction
            filename = os.path.basename(url)
            parameters = ["fastq"]
            if is_paired:
                parameters.extend(["-1", "/data/fastq/{}_1.fq".format(filename),
                                  "-2" "/data/fastq/{}_2.fq".format(filename),
                                  "-0" "/data/bam/{}_unpaired.bam".format(filename)])
            else:
                parameters.extend(["-1", "/data/fastq/{}.fq".format(filename)])
            parameters.extend(["/data/bam/{}".format(filename)])
            # extract
            job.fileStore.logToMaster("Running {} with parameters: {}".format(DOCKER_SAMTOOLS, parameters))
            dockerCall(job, tool=DOCKER_SAMTOOLS,
                       workDir=work_dir, parameters=parameters)
            # validate the unpaired reads
            unpaired_size = os.stat(os.path.join(bam_location, "{}_unpaired.bam".format(filename))).st_size
            if unpaired_size != 0:
                #todo throw usererror?
                job.fileStore.logToMaster("Got non-zero length ({}) for unmatched reads in file \"{}\""
                                          .format(unpaired_size, filename))
        # just need to download
        elif url.endswith(FQ_FORMATS) and not url.endswith(TAR_FORMATS):
            download_url(job, url=url, work_dir=fastq_location)
        # should never happen
        else:
            raise UserError('PROGRAMMER ERROR: URL "{}" with type "{}" did not match expected cases for downloading!'
                            .format(url, file_type))

    # now we need to concatinate them together
    all_fq_files = os.listdir(fastq_location)
    if len(all_fq_files) == 0:
        raise UserError("Found no fastq files after preparation for concatination.")
    output = None

    # for handling paired sample data
    if is_paired:
        #prep
        front_files = _files_with_suffix(all_fq_files, PAIRED_FRONT_SUFFIXES)
        back_files = _files_with_suffix(all_fq_files, PAIRED_BACK_SUFFIXES)
        front_files.sort()
        back_files.sort()

        #sanity check same number of files
        if len(front_files ) != len(back_files):
            raise UserError('Found inconsistent number of files ({} / {}) for "paired" input: {}'
                            .format(len(front_files), len(back_files), [os.path.basename(f) for f in all_fq_files]))
        #sanity check not missing any files
        for fq in all_fq_files:
            if fq not in front_files or fq not in back_files:
                raise UserError('Found fq file without suffix "_1"/"_2" or "R1"/"R2" for "paired" input: {}'.format(fq))
        #sanity check all files have a pair and are of the same size
        total_f_size = 0
        total_b_size = 0
        for f,b in zip(front_files, back_files):
            f_size = os.stat(f).st_size
            b_size = os.stat(b).st_size
            if f[:f.rfind(PAIRED_FRONT_SUFFIXES)] != b[:b.rfind(PAIRED_BACK_SUFFIXES)]:
                raise UserError('Found mismatched paired fq files: "{}" / "{}"'
                                .format(os.path.basename(f), os.path.basename(b)))
            if f_size != b_size:
                raise UserError('Found paired fq files of different sizes: "{}" ({}) / "{}" ({})"'
                                .format(f, f_size, b, b_size))
            total_f_size += f_size
            total_b_size += b_size

        # we know that sorted order of front_files and back_files match in name and size
        # concatinate all files
        front_outfile = os.path.join(work_dir, uuid + "_1.fq")
        back_outfile = os.path.join(work_dir, uuid + "_2.fq")
        with open(front_outfile, 'w') as outfile:
            for f in front_files:
                with open(f) as infile:
                    for line in infile:
                        outfile.write(line)
        with open(back_outfile, 'w') as outfile:
            for b in back_files:
                with open(b) as infile:
                    for line in infile:
                        outfile.write(line)

        # another sanity check
        if os.stat(front_outfile).st_size != total_f_size:
            job.fileStore.logToMaster('File size {} from concatinated front files does not match expected value {}'
                  .format(os.stat(front_outfile).st_size, total_f_size))
        if os.stat(back_outfile).st_size != total_b_size:
            job.fileStore.logToMaster('File size {} from concatinated back files does not match expected value {}'
                  .format(os.stat(back_outfile).st_size, total_b_size))

        # save to output directory
        if trevdev:
            job.fileStore.logToMaster('Moving {} to output dir: {}'.format(front_outfile, config.output_dir))
            mkdir_p(config.output_dir)
            copy_files(file_paths=[front_outfile], output_dir=config.output_dir)
            job.fileStore.logToMaster('Moving {} to output dir: {}'.format(back_outfile, config.output_dir))
            copy_files(file_paths=[back_outfile], output_dir=config.output_dir)
        else:
            front_outfile_id = job.fileStore.writeGlobalFile(front_outfile)
            back_outfile_id = job.fileStore.writeGlobalFile(back_outfile)
            output = [front_outfile_id, back_outfile_id]
    # for handling non-paired sample data
    else:
        output_file = os.path.join(work_dir, uuid + ".fq")
        total_size = 0
        with open(output_file, 'w') as outfile:
            for f in all_fq_files:
                total_size += os.stat(f).st_size
                with open(f) as infile:
                    for line in infile:
                        outfile.write(line)

        # another sanity check
        if os.stat(output_file).st_size != total_size:
            job.fileStore.logToMaster('File size {} from concatinated front files does not match expected value {}'
                                      .format(os.stat(output_file).st_size, total_size))

        # save to output directory
        if trevdev:
            job.fileStore.logToMaster('Moving {} to output dir: {}'.format(output_file, config.output_dir))
            mkdir_p(config.output_dir)
            copy_files(file_paths=[outfile], output_dir=config.output_dir)
        else:
            output_file_id = job.fileStore.writeGlobalFile(output_file)
            output = output_file

    # sanity check
    if outfile is None:
        raise UserError("PROGRAMMER ERROR: missing output file after concatination!")

    # add next job
    job.addFollowOnJobFn(perform_alignment, config, output)


    pass

def perform_alignment(job, config, fq_file):
    """
    input: one (paired potentially fastq file
    run the bwa-mem container situation
    output: bam file?
    """
    job.fileStore.logToMaster("Got to perform alignment successfully with {} input"
                              .format("single" if isinstance(fq_file, str) else "paired"))

def mark_duplicates():
    """

    """
    pass

def recalibrate_quality_scores():
    """

    """
    pass

def bin_quality_scores():
    """

    """
    pass



def generate_config():
    return textwrap.dedent("""
        # RNA-seq CGL Pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline: "toil-topmed run"
        # Just Kallisto or STAR/RSEM can be run by supplying only the inputs to those tools
        #
        # URLs can take the form: http://, ftp://, file://, s3://, gnos://
        # Local inputs follow the URL convention: file:///full/path/to/input
        # S3 URLs follow the convention: s3://bucket/directory/file.txt
        #
        # Comments (beginning with #) do not need to be removed. Optional parameters left blank are treated as false.
        ##############################################################################################################
        # Required: URL {scheme} to reference genome.  You are STRONGLY encouraged to use the GRCh38DH, 1000 Genomes Project version
        reference-genome: s3://cgl-pipeline-inputs/topmed_cgl/todo/host/that/file

        # Required: Output location of sample. Can be full path to a directory or an s3:// URL
        # Warning: S3 buckets must exist prior to upload or it will fail.
        output-dir:

        # Optional: Provide a full path to a 32-byte key used for SSE-C Encryption in Amazon
        ssec:
    """.format(scheme=[x + '://' for x in SCHEMES])[1:])


def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample to be run.
        #
        #   There are 4 tab-separated columns: filetype, paired/single, UUID, URL(s) to sample
        #
        #   UUID        This should be a unique identifier for the sample to be processed
        #   filetype    Filetype of the sample. Options: "tar", "fq", "bam"
        #   paired      Indicates whether the data is paired or single-ended. Options:  "paired" or "single"
        #   URL         A URL {scheme} pointing to the sample or a full path to a directory
        #
        #   If samples are in .tar or .tar.gz format, the files must be in the root of the tarball (not in a folder)
        #   If samples are submitted as multiple fastq files, provide comma separated URLs.
        #   If samples are submitted as multiple bam files, provide comma separated URLs.
        #   Samples must have the same extension - do not mix and match gzip and non-gzipped sample pairs.
        #
        #   Samples in "tar" format must have one of these extensions: .tar .tar.gz
        #   Files in the tarballs must have one of these extensions: .fastq .fq
        #   Paired files in tarballs must follow name conventions ending in an R1/R2 or _1/_2
        #
        #   Samples in "fq" format must have one of these extensions: fastq.gz, fastq, fq.gz, fq
        #   Paired samples in "fq" format must follow name conventions ending in an R1/R2 or _1/_2.  Ordering in manifest does not matter.
        #
        #   Samples in "bam" format will be converted to fastq with picard's SamToFastq utility.
        #   The paired/single parameter will determine whether the pipleine requests single or paired output from picard
        #
        #   When multiple samples are submitted, they will be concatinated into a single (or pair of) fastq file(s) before alignment
        #
        #   If a full path to a directory is provided for a sample, every file inside needs to be a fastq(.gz).
        #   Do not have other superfluous files / directories inside or the pipeline will complain.
        #
        #   Examples of several combinations are provided below. Lines beginning with # are ignored.
        #
        #   UUID_1  fq  single  file:///path/to/first.fq.gz,file:///path/to/second.fq.gz,file:///path/to/third.fq.gz
        #   UUID_2  fq  paired  file:///first_1.fq,file://second_2.fq.gz,file:///first_2.fq,file://second_1.fq.gz
        #   UUID_3  bam single  s3:///path/to/first.unaligned.bam,s3:///path/to/second.unaligned.bam
        #   UUID_4  bam paired  file:///path/to/unaligned.bam
        #   UUID_5  tar single  http://sample-depot.com/sample.tar
        #   UUID_6  tar paired  s3://my-bucket-name/directory/paired-sample.tar.gz
        #   UUID_7  fq  single  /full/path/to/local/directory/of/fastqs/
        #   UUID_8  bam paired  /full/path/to/local/directory/of/bams/
        #
        #   Place your samples below, one per line.
        """.format(scheme=[x + '://' for x in SCHEMES])[1:])


def generate_file(file_path, generate_func):
    """
    Checks file existance, generates file, and provides message

    :param str file_path: File location to generate file
    :param function generate_func: Function used to generate file
    """
    require(not os.path.exists(file_path), file_path + ' already exists!')
    with open(file_path, 'w') as f:
        f.write(generate_func())
    print('\t{} has been generated in the current working directory.'.format(os.path.basename(file_path)))


def main():
    """
    Computational Genomics Lab, Genomics Institute, UC Santa Cruz
    Toil TOPMed pipeline

    =======================================
    Dependencies
    Curl:       apt-get install curl
    Docker:     wget -qO- https://get.docker.com/ | sh
    Toil:       pip install toil
    Boto:       pip install boto (OPTIONAL)
    """

    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')

    # Generate subparsers
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    subparsers.add_parser('generate-manifest', help='Generates an editable manifest in the current working directory.')
    subparsers.add_parser('generate', help='Generates a config and manifest in the current working directory.')

    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the RNA-seq single cell pipeline')
    group = parser_run.add_mutually_exclusive_group()
    parser_run.add_argument('--config', default=DEFAULT_CONFIG_NAME, type=str,
                            help='Path to the (filled in) config file, generated with "generate-config". '
                                 '\nDefault value: "%(default)s"')
    group.add_argument('--manifest', default=DEFAULT_MANIFEST_NAME, type=str,
                       help='Path to the (filled in) manifest file, generated with "generate-manifest". '
                            '\nDefault value: "%(default)s"')

    # If no arguments provided, print full help menu
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # Add Toil options
    Job.Runner.addToilOptions(parser_run)
    args = parser.parse_args()

    # Parse subparsers related to generation of config and manifest
    cwd = os.getcwd()
    if args.command == 'generate-config' or args.command == 'generate':
        generate_file(os.path.join(cwd, DEFAULT_CONFIG_NAME), generate_config)
    if args.command == 'generate-manifest' or args.command == 'generate':
        generate_file(os.path.join(cwd, DEFAULT_MANIFEST_NAME), generate_manifest)

    # Pipeline execution
    elif args.command == 'run':
        # sanity check
        require(os.path.exists(args.config), '{} not found. Please run '
                                             '"toil-topmed generate-config"'.format(args.config))
        require(os.path.exists(args.manifest), '{} not found and no samples provided. Please '
                                               'run "toil-topmed generate-manifest"'.format(args.manifest))

        # get samples
        samples = parse_samples(path_to_manifest=args.manifest)

        # Parse config
        parsed_config = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).iteritems()}
        config = argparse.Namespace(**parsed_config)
        config.maxCores = int(args.maxCores) if args.maxCores else sys.maxint

        # Config sanity checks
        require(config.output_dir, 'No output location specified: {}'.format(config.output_dir))
        #todo more sanity checks
        # require(config.kallisto_index,
        #         'URLs not provided for Kallisto index, so there is nothing to do!')
        # require(urlparse(config.kallisto_index).scheme in SCHEMES,
        #         'Kallisto index in config must have the appropriate URL prefix: {}'.format(SCHEMES))
        if not config.output_dir.endswith('/'):
            config.output_dir += '/'

        # Program checks
        for program in ['curl', 'docker']:
            require(next(which(program), None), program + ' must be installed on every node.'.format(program))

        # Start the workflow
        Job.Runner.startToil(Job.wrapJobFn(map_job, prepare_input, samples, config), args)


if __name__ == '__main__':
    try:
        main()
    except UserError as e:
        print(e.message, file=sys.stderr)
        sys.exit(1)
