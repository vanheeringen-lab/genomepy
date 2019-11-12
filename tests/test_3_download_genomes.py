from tempfile import mkdtemp
import genomepy
import shutil
import pytest
import os

# # Python 2
# try:
#     FileNotFoundError
# except NameError:
#     FileNotFoundError = IOError

travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"


@pytest.fixture(scope="module", params=[".gz", ".tar.gz"])
def zipped_genomes(request):
    return request.param


def test_zipped_genomes(zipped_genomes):
    """Download a gzipped and a tar.gzipped genome"""
    genome = "AgamP4" if zipped_genomes == ".gz" else "sacCer3"
    provider = "Ensembl" if zipped_genomes == ".gz" else "UCSC"
    tmp = mkdtemp()
    genomepy.install_genome(genome, provider, genome_dir=tmp)
    shutil.rmtree(tmp)


# 2019-11-07 BDGP6, BDGP6.22 and dere_caf1 currently fails on Ensembl
# @pytest.mark.xfail()
def test_ensembl_genome():
    """Test Ensembl.

    Download Drosophila genome from Ensembl and retrieve a
    specific sequence.
    """
    tmp = mkdtemp()
    genomepy.install_genome("dere_caf1", "Ensembl", genome_dir=tmp)
    g = genomepy.Genome("dere_caf1", genome_dir=tmp)
    seq = g["scaffold_4929"][61:121]
    assert (
        str(seq).upper()
        == "ACACTTGGCTAAATTCCAAGAGAACTTACAAAAACAAACACCCCCCAGTCGAAAAAGGCC"
    )
    shutil.rmtree(tmp)


def test_ucsc_genome():
    """Test UCSC.

    Download S. cerevisiae genome from UCSC and retrieve a specific sequence."""
    tmp = mkdtemp()
    genomepy.install_genome("sacCer3", "UCSC", genome_dir=tmp)
    g = genomepy.Genome("sacCer3", genome_dir=tmp)
    seq = g["chrIV"][1337000:1337020]
    assert str(seq) == "TTTGGTTGTTCCTCTTCCTT"


# 2019-11-12 Release_6_plus_ISO1_MT currently fails on NCBI
# @pytest.mark.xfail()
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
