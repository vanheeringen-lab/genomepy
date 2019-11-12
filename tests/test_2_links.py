import genomepy
import pytest
import os

# Python 2
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"


def test_ensembl_genome_download_links():
    """Test Ensembl FTP links for various genomes

    These genomes are hosted on ftp.ensembl.org.
    """
    p = genomepy.provider.ProviderBase.create("Ensembl")

    for genome in ["GRCz11"]:  # "GRCh38.p13"
        p.get_genome_download_link(genome)


# @pytest.mark.skipif(travis, reason="FTP does not work on Travis")
def test_ensemblgenomes_genome_download_links():
    """Test Ensembl FTP links for various genomes

    These genomes are hosted on ftp.ensemblgenomes.org.
    """
    p = genomepy.provider.ProviderBase.create("Ensembl")

    for genome in ["AgamP4", "WBcel235"]:  # "R64-1-1"
        p.get_genome_download_link(genome)


@pytest.fixture(scope="module", params=["primary_assembly", "toplevel"])
def assembly(request):
    return request.param


@pytest.fixture(scope="module", params=["hard", "soft", "unmasked"])
def masking(request):
    return request.param


@pytest.fixture(scope="module", params=["54", None])
def release_version(request):
    return request.param


def test_ensembl_download_options(assembly, masking, release_version):
    """Test masking, as well as Ensembl specific options for assembly level and release versions"""
    p = genomepy.provider.ProviderBase.create("Ensembl")

    mask = masking
    toplevel = False if assembly == "primary_assembly" else True
    version = release_version

    p.get_genome_download_link(
        "GRCh38.p13", mask=mask, toplevel=toplevel, version=version
    )


def test_ucsc_genome_download_links(masking):
    """Test UCSC FTP links for various genomes

    Also test masking (unmasked should be ignored)."""
    p = genomepy.provider.ProviderBase.create("UCSC")

    for genome in ["sacCer3", "hg38"]:
        p.get_genome_download_link(genome, mask=masking)


def test_ncbi_genome_download_links(masking):
    """Test NCBI FTP links for various genomes

    Also test masking (should be ignored)."""
    p = genomepy.provider.ProviderBase.create("NCBI")

    for genome in ["Charlie1.0", "GRCh38.p13"]:
        p.get_genome_download_link(genome, mask=masking)


def test_bad_url():
    """Test URL.

    Try to download a non-genome url link.
    """
    with pytest.raises(ValueError):
        genomepy.install_genome("www.google.com", "url")


def test_nonexisting_url():
    """Test URL.

    Try to download a non-genome url link.
    """
    with pytest.raises(ValueError):
        genomepy.install_genome("this is not an url", "url")
