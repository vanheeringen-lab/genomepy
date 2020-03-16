import genomepy
import shutil
import gzip
import os
import pytest
from tempfile import mkdtemp
from platform import system

travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"
linux = system() == "Linux"


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


def test_ensembl_annotation(localname=None):
    """Test Ensembl annotation

    This annotation is hosted on https://ftp.ensembl.org.
    """
    tmp = mkdtemp()
    p = genomepy.provider.ProviderBase.create("Ensembl")

    # Only test on vertebrates as these are downloaded over HTTPS.
    # All others are downloaded over FTP, which is unreliable on Travis.
    for name, version in [("KH", 98)]:
        p.download_annotation(name, tmp, localname=localname, version=version)

        localname = genomepy.utils.get_localname(name, localname)
        gtf = os.path.join(tmp, localname, localname + ".annotation.gtf.gz")
        validate_gzipped_gtf(gtf)

        bed = os.path.join(tmp, localname, localname + ".annotation.bed.gz")
        validate_gzipped_bed(bed)

    shutil.rmtree(tmp)


@pytest.mark.skipif(travis and linux, reason="FTP does not work on Travis-Linux")
def test_ensemblgenomes_annotation(localname=None):
    """Test Ensembl annotation

    This annotation is hosted on ftp.ensemblgenomes.org.
    """
    tmp = mkdtemp()
    p = genomepy.provider.ProviderBase.create("Ensembl")

    for name, version in [("TAIR10", 45)]:
        p.download_annotation(name, tmp, localname=localname, version=version)

        localname = genomepy.utils.get_localname(name, localname)
        gtf = os.path.join(tmp, localname, localname + ".annotation.gtf.gz")
        validate_gzipped_gtf(gtf)

        bed = os.path.join(tmp, localname, localname + ".annotation.bed.gz")
        validate_gzipped_bed(bed)

    shutil.rmtree(tmp)


def test_UCSC_annotation(localname=None):
    """Test UCSC annotation"""
    tmp = mkdtemp()
    p = genomepy.provider.ProviderBase.create("UCSC")
    name = "sacCer3"

    p.download_annotation(name, tmp, localname=localname)

    localname = genomepy.utils.get_localname(name, localname)
    gtf = os.path.join(tmp, localname, localname + ".annotation.gtf.gz")
    validate_gzipped_gtf(gtf)

    bed = os.path.join(tmp, localname, localname + ".annotation.bed.gz")
    validate_gzipped_bed(bed)

    shutil.rmtree(tmp)


def test_NCBI_annotation(localname=None):
    """Test NCBI annotation"""
    tmp = mkdtemp()
    p = genomepy.provider.ProviderBase.create("NCBI")
    name = "ASM2732v1"

    p.download_annotation(name, tmp, localname=localname)

    localname = genomepy.utils.get_localname(name, localname)
    gtf = os.path.join(tmp, localname, localname + ".annotation.gtf.gz")
    validate_gzipped_gtf(gtf)

    bed = os.path.join(tmp, localname, localname + ".annotation.bed.gz")
    validate_gzipped_bed(bed)

    shutil.rmtree(tmp)


def test_url_annotation(localname=None):
    """Test url annotation"""
    tmp = mkdtemp()
    p = genomepy.provider.ProviderBase.create("URL")
    name = "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XT9_1.fa.gz"

    p.download_annotation(name, tmp, localname=localname)

    localname = genomepy.utils.get_localname(name, localname)
    gtf = os.path.join(tmp, localname, localname + ".annotation.gtf.gz")
    validate_gzipped_gtf(gtf)

    bed = os.path.join(tmp, localname, localname + ".annotation.bed.gz")
    validate_gzipped_bed(bed)

    shutil.rmtree(tmp)
