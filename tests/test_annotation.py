from tempfile import mkdtemp
import genomepy
import shutil
import gzip
import os
import pytest

# Python 2
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"


def setup():
    pass


def teardown():
    pass


def validate_gzipped_gtf(fname):
    assert os.path.exists(fname)
    with gzip.open(fname, "r") as f:
        for line in f:
            line = line.decode()
            if line.startswith("#"):
                continue
            vals = line.split("\t")
            assert 9 == len(vals)
            int(vals[3]), int(vals[4])
            break


def validate_gzipped_bed(fname):
    assert os.path.exists(fname)
    with gzip.open(fname, "r") as f:
        for line in f:
            line = line.decode()
            if line.startswith("#"):
                continue
            vals = line.split("\t")
            assert 12 == len(vals)
            int(vals[1]), int(vals[2])
            break


def test_ensembl_annotation():
    """Test Ensembl annotation

    This annotation is hosted on ftp.ensembl.org.
    """
    tmp = mkdtemp()
    p = genomepy.provider.ProviderBase.create("Ensembl")

    for name, version in [("Turkey_2.01", 98)]:
        p.download_annotation(name, tmp, localname=name, version=version)

        gtf = os.path.join(tmp, name, name + ".annotation.gtf.gz")
        validate_gzipped_gtf(gtf)

        bed = os.path.join(tmp, name, name + ".annotation.bed.gz")
        validate_gzipped_bed(bed)

    shutil.rmtree(tmp)


# @pytest.mark.skipif(travis, reason="FTP does not work on Travis")
def test_ensemblgenomes_annotation():
    """Test Ensembl annotation

    This annotation is hosted on ftp.ensemblgenomes.org.
    """
    tmp = mkdtemp()
    p = genomepy.provider.ProviderBase.create("Ensembl")

    for name, version in [("TAIR10", 45)]:
        p.download_annotation(name, tmp, version=version)

        gtf = os.path.join(tmp, name, name + ".annotation.gtf.gz")
        validate_gzipped_gtf(gtf)

        bed = os.path.join(tmp, name, name + ".annotation.bed.gz")
        validate_gzipped_bed(bed)

    shutil.rmtree(tmp)


def test_UCSC_annotation():
    """Test UCSC annotation"""
    tmp = mkdtemp()
    p = genomepy.provider.ProviderBase.create("UCSC")
    name = "sacCer3"

    p.download_annotation(name, tmp, localname=name)

    gtf = os.path.join(tmp, name, name + ".annotation.gtf.gz")
    validate_gzipped_gtf(gtf)

    bed = os.path.join(tmp, name, name + ".annotation.bed.gz")
    validate_gzipped_bed(bed)

    shutil.rmtree(tmp)


def test_NCBI_annotation():
    """Test NCBI annotation"""
    tmp = mkdtemp()
    p = genomepy.provider.ProviderBase.create("NCBI")
    name = "ASM69910v1"

    p.download_annotation(name, tmp, localname=name)

    gtf = os.path.join(tmp, name, name + ".annotation.gtf.gz")
    validate_gzipped_gtf(gtf)

    bed = os.path.join(tmp, name, name + ".annotation.bed.gz")
    validate_gzipped_bed(bed)

    shutil.rmtree(tmp)
