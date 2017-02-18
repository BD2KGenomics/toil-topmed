import logging
import os
from toil_lib import require, UserError
from toil_topmed.topmed_cgl_pipeline import parse_samples

# def parse_samples():
#     return ""

log = logging.getLogger(__name__)



def test_parse_samples_valid_cases(tmpdir):

    # fq validation
    out1 = _test_successful_parse_sample(tmpdir, "uuid", "fq", "single",
            ["file:///file.fq", "http://host/file.fastq", "s3:///bucket/file.fq.gz", "ftp://host/file.fastq.gz"])
    out2 = _test_successful_parse_sample(tmpdir, "uuid", "fq", "paired",
            ["file:///file_1.fq", "file:///file_2.fq", "file:///fileR1.fq", "file:///fileR2.fq"])

    # tar validation
    out3 = _test_successful_parse_sample(tmpdir, "uuid", "tar", "single",
            ["file:///tarball.tar"])
    out4 = _test_successful_parse_sample(tmpdir, "uuid", "tar", "paired",
            ["s3:///tarball.tar.gz"])
    out5 = _test_successful_parse_sample(tmpdir, "uuid", "tar", "single",
            ["http://host/tarball.tar"])
    out6 = _test_successful_parse_sample(tmpdir, "uuid", "tar", "paired",
            ["ftp://host/tarball.tar.gz"])

    # bam validation
    out7 = _test_successful_parse_sample(tmpdir, "uuid", "bam", "single",
            ["file:///file.bam", "http://host/file.bam", "s3:///bucket/file.bam", "ftp://host/url.bam"])
    out8 = _test_successful_parse_sample(tmpdir, "uuid", "bam", "paired",
            ["file:///file.bam", "file:///file2.bam", "file:///file3.bam"])


    # directory on filesystem
    directory = os.path.join(str(tmpdir), "test")
    os.mkdir(directory)
    filenames = ["test1.fastq.gz", "test2.fastq", "test3.fq.gz", "test4.fq"]
    [open(os.path.join(directory, f), 'a').close() for f in filenames]
    expected_filenames = ["file://" + os.path.join(directory, f) for f in filenames]
    out9 = _test_successful_parse_sample(tmpdir, "uuid", "fq", "single", [directory], expected_filenames)

    # all inputs
    input = [out1[0], out2[0], out3[0], out4[0], out5[0], out6[0], out7[0], out8[0], out9[0]]
    for i in input:
        i[3] = ",".join(i[3])
    manifest_location = _generate_manifest(tmpdir, input)
    output = parse_samples(manifest_location)
    require(len(output) == 9, "expected to have nine outputs for full parse samples test")

def test_parse_samples_error_cases(tmpdir):

    # scheme issues
    _test_failed_parse_sample(tmpdir, "uuid", "fq", "single", ["badscheme:///file.fq"],
                              "expected error for bad scheme")
    _test_failed_parse_sample(tmpdir, "uuid", "bam", "single", ["badscheme:///file.bam"],
                              "expected error for bad scheme")
    _test_failed_parse_sample(tmpdir, "uuid", "tar", "single", ["badscheme:///file.tar"],
                              "expected error for bad scheme")

    # extension issues
    _test_failed_parse_sample(tmpdir, "uuid", "fq", "single", ["file:///file.badextension"],
                              "expected error for bad extension")
    _test_failed_parse_sample(tmpdir, "uuid", "bam", "single", ["file:///file.badextension"],
                              "expected error for bad extension")
    _test_failed_parse_sample(tmpdir, "uuid", "tar", "single", ["file:///file.badextension"],
                              "expected error for bad extension")

    # paired issues
    _test_failed_parse_sample(tmpdir, "uuid", "fq", "paired", ["file:///file.fq", "file:///anotherfile.fq"],
                              "expected error for inconsistent file pairs")

    # tar issue
    _test_failed_parse_sample(tmpdir, "uuid", "tar", "single", ["file:///file.tar", "file:///anotherfile.tar"],
                              "expected error for multiple tars")

    # input issues
    _test_failed_parse_sample(tmpdir, "uuid", "badfiletype", "paired", ["file:///file.fq"],
                              "expected error for bad file type")
    _test_failed_parse_sample(tmpdir, "uuid", "fq", "badpairedtype", ["file:///file.fq"],
                              "expected error for bad paired type")
    manifest_location = _generate_manifest(tmpdir, ["uuid", "fq", "paired"])
    try:
        parse_samples(manifest_location)
    except UserError:
        pass
    else:
        require(False, "expected error for malformed manifest line")
    manifest_location = _generate_manifest(tmpdir, ["uuid", "fq", "paired", "file:///file.fq", "somethingelse"])
    try:
        parse_samples(manifest_location)
    except UserError:
        pass
    else:
        require(False, "expected error for malformed manifest line")


def _test_successful_parse_sample(tmpdir, uuid_in, type_in, paired_in, urls_in, expected_urls_out = None, validate=True):
    if expected_urls_out is None:
        expected_urls_out = urls_in
    input = [[uuid_in, type_in, paired_in, ",".join(urls_in)]]
    manifest_location = _generate_manifest(tmpdir, input)
    output = parse_samples(manifest_location)
    if not validate:
        return output

    require(len(output) == 1, "Expected to have single output")
    uuid_out, type_out, paired_out, urls_out = output[0]
    require(uuid_in == uuid_out, "Expected uuid to be '{}'. Got '{}'".format(uuid_in, uuid_out))
    require(type_in == type_out, "Expected type to be '{}'. Got '{}'".format(type_in, type_out))
    require(paired_in == paired_out, "Expected paired to be '{}'. Got '{}'".format(paired_in, paired_out))
    for uo in expected_urls_out:
        require(uo in urls_out, "Expected '{}' to be in output urls. Got '{}".format(uo, urls_out))
    for uo in urls_out:
        require(uo in expected_urls_out, "Unexpected url '{}' to be output urls. Expected '{}".format(uo, expected_urls_out))
    return output


def _test_failed_parse_sample(tmpdir, uuid_in, type_in, paired_in, urls_in, error_message):
    input = [[uuid_in, type_in, paired_in, ",".join(urls_in)]]
    manifest_location = _generate_manifest(tmpdir, input)
    try:
        parse_samples(manifest_location)
    except UserError:
        pass
    else:
        require(False, error_message)


def _generate_manifest(tmpdir, list_of_lines):
    manifest_location = os.path.join(str(tmpdir), "manifest-toil-rnaseqsc-test.tsv")
    # clear old manifest if it exists
    if os.path.isfile(manifest_location):
        os.remove(manifest_location)
    with open(manifest_location, 'w') as manifest:
        for line in list_of_lines:
            manifest.write("\t".join(line) + "\n")
    return manifest_location
