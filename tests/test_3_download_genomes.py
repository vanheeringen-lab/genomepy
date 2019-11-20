import genomepy
import shutil
import pytest
import os
from tempfile import mkdtemp
from platform import system
from urllib.request import URLError

travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"
linux = system() == "Linux"


@pytest.fixture(scope="module", params=[".gz", ".tar.gz"])
def zipped_genomes(request):
    return request.param


@pytest.mark.xfail(raises=URLError, reason="Travis sometimes cannot connect")
def test_zipped_genomes(zipped_genomes):
    """Download a gzipped and a tar.gzipped genome"""
    genome = "ASM14646v1" if zipped_genomes == ".gz" else "sacCer3"
    provider = "NCBI" if zipped_genomes == ".gz" else "UCSC"
    tmp = mkdtemp()
    genomepy.install_genome(genome, provider, genome_dir=tmp)
    shutil.rmtree(tmp)


@pytest.mark.skipif(travis, reason="Too slow for Travis")
@pytest.mark.xfail(raises=URLError, reason="Travis sometimes cannot connect")
@pytest.mark.slow
def test_ensembl_genome(genome="KH", provider="Ensembl", version=98):
    """Test Ensembl.

    Download smallest genome from Ensembl's HTTPS and retrieve a specific sequence.
    """
    tmp = mkdtemp()
    # Only test on vertebrates as these are downloaded over HTTPS.
    # All others are downloaded over FTP, which is unreliable on Travis.
    genomepy.install_genome(genome, provider, genome_dir=tmp, version=version)
    g = genomepy.Genome(genome, genome_dir=tmp)
    seq = g["1"][40:60]
    assert str(seq).upper() == "nnnnnnnnnnAACCCCTAAC".upper()
    shutil.rmtree(tmp)


@pytest.mark.xfail(raises=URLError, reason="Travis sometimes cannot connect")
def test_ucsc_genome(genome="sacCer3", provider="UCSC"):
    """Test UCSC.

    Download S. cerevisiae genome from UCSC and retrieve a specific sequence."""
    tmp = mkdtemp()
    genomepy.install_genome(genome, provider, genome_dir=tmp)
    g = genomepy.Genome(genome, genome_dir=tmp)
    seq = g["chrIV"][1337000:1337020]
    assert str(seq) == "TTTGGTTGTTCCTCTTCCTT"


@pytest.mark.xfail(raises=URLError, reason="Travis sometimes cannot connect")
def test_ncbi_genome(genome="ASM2732v1", provider="NCBI"):
    """Test NCBI.

    Download smallest genome from NCBI and retrieve a
    specific sequence.
    """
    tmp = mkdtemp()
    genomepy.install_genome(genome, provider, genome_dir=tmp)
    g = genomepy.Genome(genome, genome_dir=tmp)
    seq = g["ANONYMOUS"][80:107]
    assert str(seq).upper() == "ATACCTTCCTTAATACTGTTAAATTAT"
    shutil.rmtree(tmp)


@pytest.mark.xfail(raises=URLError, reason="Travis sometimes cannot connect")
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
