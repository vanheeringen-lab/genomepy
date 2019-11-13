import genomepy
import shutil
import pytest
import os
from tempfile import mkdtemp
from platform import system

travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"
linux = system() == "Linux"


@pytest.fixture(scope="module", params=[".gz", ".tar.gz"])
def zipped_genomes(request):
    return request.param


def test_zipped_genomes(zipped_genomes):
    """Download a gzipped and a tar.gzipped genome"""
    genome = "Release_6_plus_ISO1_MT" if zipped_genomes == ".gz" else "sacCer3"
    provider = "NCBI" if zipped_genomes == ".gz" else "UCSC"
    tmp = mkdtemp()
    genomepy.install_genome(genome, provider, genome_dir=tmp)
    shutil.rmtree(tmp)


@pytest.mark.skipif(travis, reason="Too slow for Travis")
@pytest.mark.slow
def test_ensembl_genome():
    """Test Ensembl.

    Download smallest genome from Ensembl's HTTPS and retrieve a specific sequence.
    """
    tmp = mkdtemp()
    # only test on vertebrates which are downloaded from HTTPS, as FTP is unreliable on Travis
    genomepy.install_genome("KH", "Ensembl", genome_dir=tmp, version=98)
    g = genomepy.Genome("KH", genome_dir=tmp)
    seq = g["1"][40:60]
    assert str(seq).upper() == "nnnnnnnnnnAACCCCTAAC".upper()
    shutil.rmtree(tmp)


def test_ucsc_genome():
    """Test UCSC.

    Download S. cerevisiae genome from UCSC and retrieve a specific sequence."""
    tmp = mkdtemp()
    genomepy.install_genome("sacCer3", "UCSC", genome_dir=tmp)
    g = genomepy.Genome("sacCer3", genome_dir=tmp)
    seq = g["chrIV"][1337000:1337020]
    assert str(seq) == "TTTGGTTGTTCCTCTTCCTT"


def test_ncbi_genome():
    """Test NCBI.

    Download Drosophila genome from NCBI and retrieve a
    specific sequence.
    """
    tmp = mkdtemp()
    genomepy.install_genome("Release 6 plus ISO1 MT", "NCBI", genome_dir=tmp)
    g = genomepy.Genome("Release_6_plus_ISO1_MT", genome_dir=tmp)
    seq = g["3L"][10637840:10637875]
    assert str(seq).upper() == "TTTGCAACAGCTGCCGCAGTGTGACCGTTGTACTG"
    shutil.rmtree(tmp)


def test_url_genome():
    """Test URL.

    Download S. cerevisiae genome directly from an url from UCSC and retrieve a
    specific sequence.
    """
    tmp = mkdtemp()
    genomepy.install_genome(
        "http://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/chromFa.tar.gz",
        "url",
        genome_dir=tmp,
        localname="url_test",
    )
    g = genomepy.Genome("url_test", genome_dir=tmp)
    assert str(g["chrI"][:12]).lower() == "gcctaagcctaa"
    shutil.rmtree(tmp)
