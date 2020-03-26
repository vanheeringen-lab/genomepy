import genomepy
import os
import shutil

from tempfile import NamedTemporaryFile, mkdtemp
from tests.test_04_download_annotation import validate_gzipped_bed, validate_gzipped_gtf


def test_sizes():
    infa = "tests/data/gap.fa"

    tmp = NamedTemporaryFile().name
    genomepy.utils.generate_fa_sizes(infa, tmp)

    result = open(tmp).read()

    assert result == "chr1\t28\nchr2\t45\nchr3\t15\n"


def test_gaps():
    infa = "tests/data/gap.fa"
    outbed = "tests/data/gap.bed"

    tmp = NamedTemporaryFile().name
    genomepy.utils.generate_gap_bed(infa, tmp)

    result = open(tmp).read()
    expect = open(outbed).read()

    assert result == expect


def test_sanitize_annotation(localname=None):
    """
    Test sanitizing of annotations

    Requires the genome, sizes_file and gtf to test
    """
    tmp = mkdtemp()
    name = "ASM2732v1"

    genomepy.functions.install_genome(
        name=name,
        provider="NCBI",
        genome_dir=tmp,
        localname=localname,
        annotation=True,
        force=True,
    )

    localname = genomepy.utils.get_localname(name, localname)
    gtf = os.path.join(tmp, localname, localname + ".annotation.gtf.gz")
    validate_gzipped_gtf(gtf)

    bed = os.path.join(tmp, localname, localname + ".annotation.bed.gz")
    validate_gzipped_bed(bed)

    shutil.rmtree(tmp)
