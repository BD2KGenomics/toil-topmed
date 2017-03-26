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
import math
import shutil
import glob
import time
import datetime
import logging

import yaml
from bd2k.util.files import mkdir_p
from bd2k.util.processes import which
from toil.job import Job
from toil.lib.docker import dockerCall, _fixPermissions
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

# resource estimatin
PREPARE_FS_TO_TO_DSK_REQ = 2.5
ALIGN_FS_TO_DSK_REQ = 6
MERGE_FS_TO_TO_DSK_REQ = 2.5
PICARD_FS_TO_TO_DSK_REQ = 2.5
GATK_FS_TO_TO_DSK_REQ = 2.5
CRAM_FS_TO_TO_DSK_REQ = 2.5
REFERENCE_SIZE = 8 * 1024 * 1024 * 1024 #8G
MAX_MEMORY = 32 * 1024 * 1024 * 1024 #32GB

# filenames
DEFAULT_CONFIG_NAME = 'config-toil-topmed.yaml'
DEFAULT_MANIFEST_NAME = 'manifest-toil-topmed.tsv'

# docker images
#todo unlatest this
DOCKER_SAMTOOLS = "quay.io/ucsc_cgl/samtools:1.3--256539928ea162949d8a65ca5c79a72ef557ce7c"
DOCKER_BWAKIT = "quay.io/ucsc_cgl/bwakit:0.7.15--ed1aeaaebf4d88ba51042da83f02ef8373a822b9"
DOCKER_GATK = "quay.io/ucsc_cgl/gatk:3.7--e931e1ca12f6d930b755e5ac9c0c4ca266370b7b"
DOCKER_PICARD = "quay.io/ucsc_cgl/picardtools:2.9.0--4d726c4a1386d4252a0fc72d49b1d3f5b50b1e23"




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
                    # valid_types = TAR_FORMATS
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

    # job prep
    config = argparse.Namespace(**vars(config))
    uuid, file_type, paired, urls = sample
    config.uuid = uuid
    config.file_type = file_type
    config.is_paired = paired
    work_dir = job.fileStore.getLocalTempDir()

    # global resource estimation
    config.cores = min(config.maxCores, multiprocessing.cpu_count())
    config.base_sample_size = config.input_file_size * len(urls)
    #todo remove this, should be specified
    config.disk = int(min(config.maxDisk, config.base_sample_size * ALIGN_FS_TO_DSK_REQ))


    # need to hand each url off to appropriate processing
    bam_ids = []
    temporary_file_ids = [] #this is to track jobstore files which can be removed
    for url in urls:
        # alignment resource estimation
        aln_mem = str(config.alignment_memory if config.alignment_memory is not None else config.memory)
        aln_cpu = config.alignment_cores if config.alignment_cores is not None else config.cores

        if file_type == 'bam':
            aln_disk = int(config.input_file_size * ALIGN_FS_TO_DSK_REQ + REFERENCE_SIZE)
            job.fileStore.logToMaster("Spawning job to handle {}".format(url))
            child_job = job.addChildJobFn(perform_alignment, config, "bam", sample_url=url, memory=aln_mem, cores=aln_cpu, disk=aln_disk)
            bam_ids.append(child_job.rv(0))
            temporary_file_ids.append(child_job.rv(1))
        elif file_type == 'tar':
            # name prep
            tar_name =  os.path.basename(url)
            tar_base = tar_name.split(".")[0]
            tar_work_dir = os.path.join(work_dir, tar_base)
            mkdir_p(tar_work_dir)

            # get and untar
            download_url(job, url=url, work_dir=tar_work_dir)
            tar_file_location = os.path.join(tar_work_dir, tar_name)
            subprocess.check_call(['tar', '-xvf', tar_file_location, '-C', tar_work_dir])

            #cleanup and post-process
            os.remove(tar_file_location)
            if os.path.isdir(os.path.join(tar_work_dir, tar_base)):
                for file in glob.glob(os.path.join(tar_work_dir, tar_base, "*")):
                    dest = os.path.join(tar_work_dir, os.path.basename(file))
                    shutil.move(file, dest)
                os.rmdir(os.path.join(tar_work_dir, tar_base))

            # handle file by type
            file_count = 0
            for bam_file in glob.glob(os.path.join(tar_work_dir, "*.bam")):
                bam_file_id = job.fileStore.writeGlobalFile(bam_file)
                job.fileStore.logToMaster("Spawning job to handle bam '{}' in '{}'".format(os.path.basename(bam_file), url))
                aln_disk = int(os.stat(bam_file).st_size * ALIGN_FS_TO_DSK_REQ + REFERENCE_SIZE)
                child_job = job.addChildJobFn(perform_alignment, config, "bam", sample_id=bam_file_id, memory=aln_mem,
                                              cores=aln_cpu, disk=aln_disk)
                bam_ids.append(child_job.rv(0))
                temporary_file_ids.append(child_job.rv(1))
                file_count += 1
            for fq_glob in [os.path.join(tar_work_dir, ("*" + fq)) for fq in FQ_FORMATS]:
                for fq in glob.glob(fq_glob):
                    raise UserError("FASTQ files in 'tar' not supported: {}".format(fq))
                    # fq_file_id = job.fileStore.writeGlobalFile(fq)
                    # child_job = job.addChildJobFn(perform_alignment, config, "fq", sample_id=fq_file_id, memory=aln_mem,
                    #                               cores=aln_cpu, disk=aln_disk)
                    # bam_ids.append(child_job.rv(0))
                    # temporary_file_ids.append(child_job.rv(1))

            #sanity check
            job.fileStore.logToMaster("Spawned {} jobs to handle {}".format(file_count, url))
            if file_count == 0:
                raise UserError("Got no files in tar: '{}'".format(url))


        # just need to download
        elif file_type == 'fq':
            # download_url(job, url=url, work_dir=fastq_location)
            #todo: how to specify read group?
            raise UserError("FQ not implemented")
        # should never happen
        else:
            raise UserError('PROGRAMMER ERROR: URL "{}" with type "{}" did not match expected cases for preprocessing'
                            .format(url, file_type))

    job.addFollowOnJobFn(merge_sam_files, config, bam_ids, temporary_file_ids,
                         memory=config.memory, cores=config.cores, disk=int(config.base_sample_size * MERGE_FS_TO_TO_DSK_REQ))


def perform_alignment(job, config, input_type, sample_url=None, sample_id=None):

    # validate input
    if sample_url is None and sample_id is None:
        raise UserError("perform_alignment invoked with neither sample_url nor sample_id")
    if sample_url is not None and sample_id is not None:
        raise UserError("perform_alignment invoked with both sample_url and sample_id")
    if input_type != "bam":
        raise UserError("input format '{}' not supported for perform_alignment".format(input_type))

    # prep
    job.fileStore.logToMaster("START_TIME:{}:{}".format(config.uuid, datetime.datetime.now()))
    start = time.time()
    work_dir = job.fileStore.getLocalTempDir()
    is_paired = config.is_paired
    is_global_file = sample_id is not None

    # download bam
    filename = "{}.in.bam".format(config.uuid)
    output_filename="{}.out.bam".format(config.uuid)
    if is_global_file:
        job.fileStore.readGlobalFile(sample_id, os.path.join(work_dir, filename))
    else:
        download_url(job, url=sample_url, work_dir=work_dir)
        filename = os.path.basename(sample_url)
        output_filename=config.uuid+"."+filename
    if not os.path.isfile(os.path.join(work_dir, filename)):
        raise UserError("Failed to download {}".format(sample_url))

    job.fileStore.logToMaster("Performing alignment on {} with {} input"
                              .format(filename, "single" if not is_paired else "paired"))

    # get fq files from bam
    fq_files, read_group_header = extract_fastq_files_from_bam(job, work_dir, filename, is_paired)

    # get reference files
    reference_location = os.path.join(work_dir, "reference")
    mkdir_p(reference_location)
    reference_name = download_reference(job, reference_location, config.reference)

    # bwa
    cores = config.alignment_cores if config.alignment_cores is not None else config.cores
    bwa_cmd = 'bwa mem -t {} -K 100000000 -R {} -Y {} {} '.format(
        cores, read_group_header, os.path.join("/data/reference", os.path.basename(reference_name)),
        os.path.join("/data", os.path.basename(fq_files[0])))
    if len(fq_files) > 1:
        bwa_cmd += os.path.join("/data", os.path.basename(fq_files[1]))
    # other commands
    samblaster_cmd = "samblaster --addMateTags"
    samtools_sam2bam_cmd = "samtools view -Sb -"
    #todo set cores on sort cmd.  split cores?
    samtools_sort_cmd = "samtools sort -@ {} -T {} -o {}".format(cores, config.uuid + ".sorting",
                                                                 os.path.join("/data", output_filename))

    # modification for the docker file I use
    bwa_cmd = "/opt/bwa.kit/"+bwa_cmd
    samblaster_cmd = "/opt/samblaster-v.0.1.24/"+samblaster_cmd
    samtools_sam2bam_cmd = "/opt/bwa.kit/"+samtools_sam2bam_cmd
    samtools_sort_cmd = "/opt/bwa.kit/"+samtools_sort_cmd

    # run the command
    params = [bwa_cmd.split(" "), samblaster_cmd.split(), samtools_sam2bam_cmd.split(), samtools_sort_cmd.split()]
    job.fileStore.logToMaster("Running {} with parameters: {}".format(DOCKER_BWAKIT, params))
    dockerCall(job, tool=DOCKER_BWAKIT, workDir=work_dir, parameters=params)

    # sanity check
    output_bam_location = os.path.join(work_dir, output_filename)
    if not os.path.isfile(output_bam_location):
        raise UserError('Aligned BAM file "{}" does not exist'.format(output_bam_location))

    _fixPermissions(DOCKER_BWAKIT, work_dir) #todo this is a temporary fix

    # save to output directory (for development)
    if config.save_intermediate_files:
        job.fileStore.logToMaster('Moving {} to output dir: {}'.format(output_bam_location, config.output_dir))
        copy_files(file_paths=[output_bam_location], output_dir=config.output_dir)

    # it worked
    aligned_bam_id = job.fileStore.writeGlobalFile(output_bam_location)
    job.fileStore.logToMaster("TIME:{}:perform_alignment:{}".format(config.uuid, time.time() - start))

    # return for next job (return sample id even if null)
    return aligned_bam_id, sample_id


def extract_fastq_files_from_bam(job, work_dir, bam_name, is_paired, get_readgroup_header=True):

    # get readgroup header
    readgroup_header = None
    if get_readgroup_header:
        headers_name = bam_name + ".headers.txt"
        headers_location = os.path.join(work_dir, headers_name)
        cmd = ["view", "-H", os.path.join("/data", bam_name)]
        with open(headers_location, 'w') as file:
            dockerCall(job, tool=DOCKER_SAMTOOLS, workDir=work_dir, parameters=cmd, outfile=file)

        # read file
        if not os.path.isfile(headers_location):
            raise UserError("Found no headers for {}".format(bam_name))
        with open(headers_location, 'r') as f:
            for l in f:
                if l.startswith("@RG"):
                    # check for multiple @RG headers
                    if readgroup_header is not None and l != readgroup_header:
                        job.fileStore.logToMaster("Found extra @RG header for file \"{}\": \"{}\""
                                        .format(bam_name, "\\t".join(l.rstrip().split("\t"))), level=logging.WARN)
                    else:
                        readgroup_header = l.rstrip()
                        # fix readgroup #todo maybe move this to alignment?
                        readgroup_header = "\\t".join(readgroup_header.split("\t"))
                        job.fileStore.logToMaster("Found @RG header for file \"{}\": \"{}\""
                                                  .format(bam_name, readgroup_header))
        if readgroup_header is None:
            raise UserError("Found no @RG header for file \"{}\"".format(bam_name))


    # prep for extraction into FQ
    fastq_files=[]
    bam_base = os.path.splitext(bam_name)[0]
    docker_params = ["fastq"]
    if is_paired:
        fastq_files.append("{}_1.fq".format(bam_base))
        fastq_files.append("{}_2.fq".format(bam_base))
        fastq_files.append("{}_unpaired.bam".format(bam_base))
        docker_params.extend(["-1", "{}".format(fastq_files[0]),
                              "-2", os.path.join("/data", fastq_files[1]),
                              "-0", os.path.join("/data", fastq_files[2])])
    else:
        fastq_files.append("{}.fq".format(bam_base))
        docker_params.extend(["-1", os.path.join("/data", fastq_files[0])])
    docker_params.extend([os.path.join("/data", bam_name)])
    # extract
    job.fileStore.logToMaster("Running {} with parameters: {}".format(DOCKER_SAMTOOLS, docker_params))
    dockerCall(job, tool=DOCKER_SAMTOOLS,
               workDir=work_dir, parameters=docker_params)
    _fixPermissions(DOCKER_SAMTOOLS, work_dir)  # todo tmpfix

    # validate the unpaired reads
    if is_paired:
        unpaired_size = os.stat(os.path.join(work_dir, fastq_files[2])).st_size
        if unpaired_size != 0:
            #maybe this should raise an exception
            job.fileStore.logToMaster("Got non-zero length ({}) for unmatched reads in file \"{}\""
                                      .format(unpaired_size, bam_name))

    # return fastqs
    return fastq_files, readgroup_header


def merge_sam_files(job, config, aligned_bam_ids, removable_file_ids=None):
    start = time.time()

    # cleanup of unnecessary files from previous steps
    remove_intermediate_jobstore_files(job, removable_file_ids)

    #setup
    output_file_name = config.uuid + ".merged.bam"
    job.fileStore.logToMaster("Merging {} BAM files".format(len(aligned_bam_ids)))
    work_dir = job.fileStore.getLocalTempDir()

    # read file
    bam_dir = os.path.join(work_dir, "bam")
    docker_bam_names = []
    os.mkdir(bam_dir)
    digits = int(math.log10(len(aligned_bam_ids))) + 1
    idx = 0
    for bam in aligned_bam_ids:
        bam_name = config.uuid + '.{0:0{digits}}.bam'.format(idx, digits=digits)
        job.fileStore.readGlobalFile(bam, os.path.join(bam_dir, bam_name))
        docker_bam_names.append(os.path.join("/data", "bam", bam_name))
        idx += 1

    # Call docker image
    # -c/-p help conflicting RG or PG headers
    params = ["merge", "-c", "--threads", str(config.cores), "-p", os.path.join("/data", output_file_name)]
    params.extend(docker_bam_names)
    job.fileStore.logToMaster("Calling {} with params: {}".format(DOCKER_SAMTOOLS, params))
    dockerCall(job, tool=DOCKER_SAMTOOLS,
               workDir=work_dir, parameters=params)
    _fixPermissions(DOCKER_SAMTOOLS, work_dir)  # todo tmpfix

    # sanity check
    output_file_location = os.path.join(work_dir, output_file_name)
    if not os.path.isfile(output_file_location):
        raise UserError('Merged BAM file "{}" does not exist'.format(output_file_location))

    # save to output directory (for development)
    if config.save_intermediate_files:
        job.fileStore.logToMaster('Moving {} to output dir: {}'.format(output_file_location, config.output_dir))
        copy_files(file_paths=[output_file_location], output_dir=config.output_dir)

    # save output
    output_id = job.fileStore.writeGlobalFile(output_file_location)

    # add next job
    job.addFollowOnJobFn(mark_duplicates, config, output_id, aligned_bam_ids,
                         memory=config.memory, cores=config.cores,
                         disk=int(config.base_sample_size * PICARD_FS_TO_TO_DSK_REQ))
    job.fileStore.logToMaster("TIME:{}:merge_sam_files:{}".format(config.uuid, time.time() - start))


def _get_default_docker_params(work_dir):
    return ['--rm','--log-driver','none', '-v', '{}:/data'.format(work_dir)]

def mark_duplicates(job, config, aligned_bam_id, removable_file_ids=None):
    start = time.time()

    # cleanup of unnecessary files from previous steps
    remove_intermediate_jobstore_files(job, removable_file_ids)

    #prep
    duped_bam_name = config.uuid + ".duped.bam"
    deduped_bam_name = config.uuid + ".deduped.bam"
    job.fileStore.logToMaster('Marking duplicates on {}'.format(duped_bam_name))
    work_dir = job.fileStore.getLocalTempDir()
    mkdir_p(os.path.join(work_dir, ".java_tmp"))

    # read file
    job.fileStore.readGlobalFile(aligned_bam_id, os.path.join(work_dir, duped_bam_name))

    # Call docker image
    #todo -Djava.io.tmpdir=/tmp
    # picard MarkDuplicates I={} O={} M={} ASSUME_SORT_ORDER=coordinate
    metrics_filename = config.uuid + ".deduplication_metrics.txt"
    params = ["MarkDuplicates",
              "I={}".format(os.path.join("/data", duped_bam_name)),
              "O={}".format(os.path.join("/data", deduped_bam_name)),
              "M={}".format(os.path.join("/data", metrics_filename)),
              "ASSUME_SORT_ORDER=coordinate"]
    job.fileStore.logToMaster("Calling {} with params: {}".format(DOCKER_PICARD, params))

    docker_params = _get_default_docker_params(work_dir)
    docker_params.extend(['-e', 'JAVA_OPTS=-Djava.io.tmpdir=/data/.java_tmp'])
    dockerCall(job, tool=DOCKER_PICARD, workDir=work_dir, parameters=params, dockerParameters=docker_params)

    # verify output
    deduped_bam_location = os.path.join(work_dir, deduped_bam_name)
    metrics_location = os.path.join(work_dir, metrics_filename)
    save_metrics = False
    if not os.path.isfile(metrics_location):
        job.fileStore.logToMaster("Deduplication metrics file was not found: {}".format(metrics_location))
    elif os.stat(os.path.join(metrics_location)).st_size != 0:
        #todo user error?  just save it?
        job.fileStore.logToMaster("Deduplication metrics file has nonzero size: {}".format(metrics_location))
        save_metrics = True

    # sanity check
    if not os.path.isfile(deduped_bam_location):
        raise UserError('Deduplicated BAM file "{}" does not exist'.format(deduped_bam_location))
    deduped_bam_id = job.fileStore.writeGlobalFile(deduped_bam_location)

    # save to output directory (for development)
    if config.save_intermediate_files:
        outfiles = [deduped_bam_location]
        if save_metrics:
            outfiles.append(metrics_location)
        job.fileStore.logToMaster('Moving {} to output dir: {}'.format(outfiles, config.output_dir))
        copy_files(file_paths=outfiles, output_dir=config.output_dir)

    # add next job
    job.addFollowOnJobFn(recalibrate_quality_scores, config, deduped_bam_id, [aligned_bam_id],
                         memory=config.memory, cores=config.cores,
                         disk=int(config.base_sample_size * GATK_FS_TO_TO_DSK_REQ))
    job.fileStore.logToMaster("TIME:{}:mark_duplicates:{}".format(config.uuid, time.time() - start))


def recalibrate_quality_scores(job, config, input_bam_id, removable_file_ids=None):
    start = time.time()

    # cleanup of unnecessary files from previous steps
    remove_intermediate_jobstore_files(job, removable_file_ids)

    #prep
    recalibration_bam_name = config.uuid + ".recalibration.bam"
    recalibration_report_name = config.uuid + ".recalibration.grp"
    job.fileStore.logToMaster('Recalibrating quality scores on {}'.format(recalibration_bam_name))

    # read file
    work_dir = job.fileStore.getLocalTempDir()
    job.fileStore.readGlobalFile(input_bam_id, os.path.join(work_dir, recalibration_bam_name))

    # get reference
    reference_dir = os.path.join(work_dir, "reference")
    mkdir_p(reference_dir)
    reference_location = download_reference(job, reference_dir, config.reference)

    # get variants
    variant_dir = os.path.join(work_dir, "variants")
    mkdir_p(variant_dir)
    variants = download_variant(job, variant_dir, config.variants)

    # get bam index
    bai_location = index_bam(job, work_dir, recalibration_bam_name)

    # build params for docker call
    params = ["-T", "BaseRecalibrator",
              "-R", os.path.join("/data/reference", os.path.basename(reference_location)),
              "-I", os.path.join("/data", recalibration_bam_name),
              "-o", os.path.join("/data", recalibration_report_name)]
    for variant in variants:
        params.extend(["-knownSites", os.path.join("/data/variants", variant)])
    optional_params = ["--downsample_to_fraction", ".1", "-L", "chr1", "-L", "chr2", "-L", "chr3", "-L", "chr4",
                       "-L", "chr5", "-L", "chr6", "-L", "chr7", "-L", "chr8", "-L", "chr9", "-L", "chr10",
                       "-L", "chr11", "-L", "chr12", "-L", "chr13", "-L", "chr14", "-L", "chr15",
                       "-L", "chr16", "-L", "chr17", "-L", "chr18", "-L", "chr19", "-L", "chr20",
                       "-L", "chr21", "-L", "chr22", "--interval_padding", "200", "-rf", "BadCigar",
                       "-nct", str(config.cores),
                       "--preserve_qscores_less_than", "6",
                       "--disable_auto_index_creation_and_locking_when_reading_rods", "--disable_bam_indexing",
                       "--useOriginalQualities"]
    params.extend(optional_params)

    # call docker: type(parameters) == list of lists => execute command as shell
    job.fileStore.logToMaster("Calling {} with command: {}".format(DOCKER_GATK, params))
    dockerCall(job, tool=DOCKER_GATK, workDir=work_dir, parameters=params)

    # sanity check
    recalibration_report_location = os.path.join(work_dir, recalibration_report_name)
    if not os.path.isfile(recalibration_report_location):
        raise UserError('Recalibration report file "{}" does not exist'.format(recalibration_report_location))
    if os.stat(recalibration_report_location).st_size == 0:
        raise UserError('Recalibration report file "{}" is empty'.format(recalibration_report_location))
    subprocess.check_call(["chmod", "+w", recalibration_report_location])

    # save file for next step
    recalibration_report_id = job.fileStore.writeGlobalFile(recalibration_report_location)
    bam_index_id = job.fileStore.writeGlobalFile(bai_location)

    # save to output directory (for development)
    if config.save_intermediate_files:
        job.fileStore.logToMaster('Moving {} to output dir: {}'
                                  .format([recalibration_report_location, bai_location], config.output_dir))
        copy_files(file_paths=[recalibration_report_location], output_dir=config.output_dir)

    # add next job
    job.addFollowOnJobFn(bin_quality_scores, config, input_bam_id, bam_index_id, recalibration_report_id,
                         memory=config.memory, cores=config.cores,
                         disk=int(config.base_sample_size * GATK_FS_TO_TO_DSK_REQ))
    job.fileStore.logToMaster("TIME:{}:recalibrate_quality_scores:{}".format(config.uuid, time.time() - start))

#todo: merge these

def bin_quality_scores(job, config, input_bam_id, bam_index_id, bsqr_report_id, removable_file_ids=None):
    start = time.time()

    # cleanup of unnecessary files from previous steps
    remove_intermediate_jobstore_files(job, removable_file_ids)

    #prep
    input_bam_name = config.uuid + ".fullqs.bam"
    input_bam_index_name = input_bam_name + ".bai"
    bqsr_report_name = config.uuid + ".bqsr.txt"
    output_bam_name = config.uuid + ".binnedqs.bam"
    job.fileStore.logToMaster('Binning quality scores on {}'.format(input_bam_name))

    # read file
    work_dir = job.fileStore.getLocalTempDir()
    input_bam_location = os.path.join(work_dir, input_bam_name)
    input_bam_index_location = os.path.join(work_dir, input_bam_index_name)
    bsqr_report_location = os.path.join(work_dir, bqsr_report_name)
    job.fileStore.readGlobalFile(input_bam_id, input_bam_location)
    job.fileStore.readGlobalFile(bam_index_id, input_bam_index_location)
    job.fileStore.readGlobalFile(bsqr_report_id, bsqr_report_location)

    # get reference
    reference_dir = os.path.join(work_dir, "reference")
    mkdir_p(reference_dir)
    reference_name = download_reference(job, reference_dir, config.reference)

    #prep commands
    params = ["-T", "PrintReads",
              "-R", os.path.join("/data/reference", os.path.basename(reference_name)),
              "-I", os.path.join("/data", input_bam_name),
              "-o", os.path.join("/data", output_bam_name),
              "-BQSR", os.path.join("/data", bqsr_report_name),
              "-SQQ", "10", "-SQQ", "20", "-SQQ", "30"]
    optional_params = ["--globalQScorePrior", "-1.0",
                       "--preserve_qscores_less_than", "6",
                       "--disable_indel_quals",
                       "--useOriginalQualities",
                       "-rf", "BadCigar",
                       "-nct", str(config.cores),
                       "--emit_original_quals"]
                       # "--createOutputBamMD5",
                       # "--addOutputSAMProgramRecord",
    params.extend(optional_params)

    # call out to docker
    output_bam_location = os.path.join(work_dir, output_bam_name)
    job.fileStore.logToMaster("Calling {} with params: {}".format(DOCKER_GATK, params))
    dockerCall(job, tool=DOCKER_GATK, workDir=work_dir, parameters=params)

    # sanity check
    if not os.path.isfile(output_bam_location):
        raise UserError('Binned-quality-scores BAM file "{}" does not exist'.format(output_bam_location))
    output_bam_id = job.fileStore.writeGlobalFile(output_bam_location)

    # save to output directory (for development)
    if config.save_intermediate_files:
        job.fileStore.logToMaster('Moving {} to output dir: {}'.format(output_bam_location, config.output_dir))
        copy_files(file_paths=[output_bam_location], output_dir=config.output_dir)

    # add next job
    job.addFollowOnJobFn(convert_to_cram_and_validate, config, output_bam_id, [input_bam_id, bam_index_id, bsqr_report_id],
                         memory=config.memory, cores=config.cores,
                         disk=int(config.base_sample_size * CRAM_FS_TO_TO_DSK_REQ))
    job.fileStore.logToMaster("TIME:{}:bin_quality_scores:{}".format(config.uuid, time.time() - start))


def convert_to_cram_and_validate(job, config, input_bam_id, removable_file_ids=None):
    start = time.time()

    # cleanup of unnecessary files from previous steps
    remove_intermediate_jobstore_files(job, removable_file_ids)

    #prep
    input_bam_name = config.uuid + ".bam"
    output_cram_name = config.uuid + ".cram"
    validation_file_name = config.uuid + ".validation.log"
    job.fileStore.logToMaster('Converting BAM to CRAM and validating: {}'.format(input_bam_name))

    # read file
    work_dir = job.fileStore.getLocalTempDir()
    job.fileStore.readGlobalFile(input_bam_id, os.path.join(work_dir, input_bam_name))

    # get reference
    reference_dir = os.path.join(work_dir, "reference")
    mkdir_p(reference_dir)
    reference_name = download_reference(job, reference_dir, config.reference)

    # convert to cram
    params = ['view', '-C',
              '-T', os.path.join("/data/reference", reference_name),
              '-o', os.path.join("/data", output_cram_name),
              os.path.join("/data", input_bam_name)]
    job.fileStore.logToMaster("Calling {} with params: {}".format(DOCKER_SAMTOOLS, params))
    dockerCall(job, tool=DOCKER_SAMTOOLS, workDir=work_dir, parameters=params)

    # sanity check
    output_cram_location = os.path.join(work_dir, output_cram_name)
    # sanity check
    if not os.path.isfile(output_cram_location):
        raise UserError('Output CRAM file "{}" does not exist'.format(output_cram_location))

    #todo validate

    # save to output
    output_bam_location = os.path.join(work_dir, input_bam_name)
    output_files = [output_cram_location, output_bam_location]
    if urlparse(config.output_dir).scheme == 's3':
        for output_file in output_files:
            job.fileStore.logToMaster('Uploading {} to S3: {}'.format(config.uuid, config.output_dir))
            s3am_upload(fpath=output_file, s3_dir=config.output_dir, num_cores=config.cores)
    else:
        job.fileStore.logToMaster('Moving {} to output dir: {}'.format(config.uuid, config.output_dir))
        mkdir_p(config.output_dir)
        copy_files(file_paths=[os.path.join(work_dir, config.uuid + '.tar.gz')], output_dir=config.output_dir)
    job.fileStore.logToMaster("TIME:{}:validate_output:{}".format(config.uuid, time.time() - start))
    job.fileStore.logToMaster("END_TIME:{}:{}".format(config.uuid, datetime.datetime.now()))


def remove_intermediate_jobstore_files(job, file_id_list):
    # for cases where there is nothing to remove: let this function handle it
    if file_id_list is None or len(file_id_list) == 0: return
    # one case has file_id_list as list of list of files
    removed_files_count = 0
    for file_id in file_id_list:
        if file_id is None:
            continue
        elif isinstance(file_id, list):
            for fid in file_id:
                job.fileStore.deleteGlobalFile(fid)
                removed_files_count += 1
        else:
            job.fileStore.deleteGlobalFile(file_id)
            removed_files_count += 1
    # log it
    job.fileStore.logToMaster("Removed {} intermediate files".format(len(file_id_list)))


def index_bam(job, work_dir, bam_name):
    start = time.time()
    bam_location = os.path.join(work_dir, bam_name)
    bai_location = "{}.bai".format(bam_location)
    job.fileStore.logToMaster("Indexing {}".format(bam_location))
    params = ["index", '-b', os.path.join("/data", bam_name)]
    dockerCall(job, DOCKER_SAMTOOLS, workDir=work_dir, parameters=params)
    if not os.path.isfile(bai_location):
        raise UserError("File not found after indexing BAM: {}".format(bai_location))
    else:
        job.fileStore.logToMaster("Index created: {}".format(bai_location))
    job.fileStore.logToMaster("TIME::index_bam:{}".format(time.time() - start))
    return bai_location


def download_reference(job, work_dir, reference_location):

    start = time.time()

    # get files and names
    download_url(job, url=reference_location, work_dir=work_dir)
    reference_tar_name = os.path.basename(reference_location)
    reference_tarfile_location = os.path.join(work_dir, reference_tar_name)

    # untar
    subprocess.check_call(['tar', '-xvf', reference_tarfile_location, '-C', work_dir])
    os.remove(reference_tarfile_location)

    # potentially move out of inner directory
    reference_base_name = reference_tar_name.split(".")[0]
    if os.path.isdir(os.path.join(work_dir, reference_base_name)):
        for file in glob.glob(os.path.join(work_dir, reference_base_name, "*")):
            dest = os.path.join(work_dir, os.path.basename(file))
            shutil.move(file, dest)
        os.rmdir(os.path.join(work_dir, reference_base_name))

    # get fq file of reference (and sanity check)
    reference_fa_name = reference_base_name + ".fa"
    if not os.path.isfile(os.path.join(work_dir, reference_fa_name)):
        raise UserError("Reference tar '{}' not in expected format: {} => [{}/]{}"
                        .format(reference_location, reference_tar_name, reference_base_name, reference_fa_name))

    job.fileStore.logToMaster("TIME::download_reference:{}".format(time.time() - start))
    return reference_fa_name


def download_variant(job, work_dir, variant_location):
    start = time.time()

    # we need this for file detection
    if len(glob.glob(os.path.join(work_dir, "*"))) != 0:
        raise UserError("download_variant must have an empty work dir")

    # get the variant tarball
    download_url(job, url=variant_location, work_dir=work_dir)
    variant_tar_name = os.path.basename(variant_location)
    variant_tarfile_location = os.path.join(work_dir, variant_tar_name)

    # untar
    subprocess.check_call(['tar', '-xvf', variant_tarfile_location, '-C', work_dir])
    os.remove(variant_tarfile_location)

    # potentially move out of inner directory
    variant_base_name = variant_tar_name.split(".")[0]
    if os.path.isdir(os.path.join(work_dir, variant_base_name)):
        for file in glob.glob(os.path.join(work_dir, variant_base_name, "*")):
            dest = os.path.join(work_dir, os.path.basename(file))
            shutil.move(file, dest)
        os.rmdir(os.path.join(work_dir, variant_base_name))

    # variants to return are prefix of another file
    variant_names = glob.glob(os.path.join(work_dir, "*"))
    variant_bases = set()
    for variant_name in variant_names:
        for vn in variant_names:
            if vn != variant_name and vn.startswith(variant_name):
                variant_bases.add(os.path.basename(variant_name))
    if len(variant_bases) == 0:
        raise UserError("Variant tar '{}' not in expected format: {} => [{}/](example.vcf[.gz], example.vcf[.gz].[tbi|idx])+"
                        .format(variant_location, variant_base_name))

    job.fileStore.logToMaster("TIME::download_variant:{}".format(time.time() - start))
    return variant_bases


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
        # Required: URL {scheme} to reference genome tarball.  Defaults to the GRCh38DH, 1000 Genomes Project version
        reference: s3://cgl-pipeline-inputs/topmed/GRCh38_full_analysis_set_plus_decoy_hla.tar.gz

        # Required: URL {scheme} to variant files tarball.  All .vcf or .vcf.gz files will be used during base recalibration
        variants: s3://cgl-pipeline-inputs/topmed/GATK_hg38_variants.tar.gz

        # Required: Output location of sample. Can be full path to a directory or an s3:// URL
        # Warning: S3 buckets must exist prior to upload or it will fail.
        output-dir:

        # Optional: Input file size in bytes (for resource estimation)
        input-file-size:

        ### Alignment Resource Management
        # Optional: CPU cores to use during alignment. default: no special alignment configuration. ex: 8
        alignment-cores:
        # Optional: Memory string to use during alignment. default: no special alignment configuration. ex: 8G
        alignment-memory:

        # Optional: Save intermediate files to the output directory (for debugging)
        save-intermediate-files: False

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
        config.maxDisk = int(args.maxDisk) if args.maxDisk else sys.maxint
        config.maxMemory = args.maxMemory if args.maxMemory else str(sys.maxint)
        config.memory = args.maxMemory if args.maxMemory else str(MAX_MEMORY)

        # Config sanity checks
        require(config.output_dir, 'No output location specified')
        if urlparse(config.output_dir).scheme != "s3":
            mkdir_p(config.output_dir)
        require(config.reference, 'No reference tarball specified')
        require(config.variants, 'No variants tarball specified')
        require(config.input_file_size, "Configuration parameter input-file-size is required")
        if config.save_intermediate_files is None or not isinstance(config.save_intermediate_files, bool):
            config.save_intermediate_files = False
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
