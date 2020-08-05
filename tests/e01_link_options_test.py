import genomepy
import pytest
import os

from platform import system
from time import sleep

travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"
linux = system() == "Linux"

skip = False
if not skip:

    @pytest.fixture(scope="module", params=["primary_assembly", "toplevel"])
    def assembly(request):
        return request.param

    @pytest.fixture(scope="module", params=["98", None])
    def release_version(request):
        return request.param

    @pytest.fixture(scope="module", params=["hard", "soft", "unmasked"])
    def masking(request):
        return request.param

    def test_ensembl_genome_download_links(assembly, masking, release_version):
        """Test Ensembl links with various options

        These genomes are hosted on ftp.ensembl.org

        Vertebrates are downloaded from HTTP.
        """
        p = genomepy.provider.ProviderBase.create("Ensembl")

        mask = masking if masking != "unmasked" else "none"
        toplevel = False if assembly == "primary_assembly" else True
        version = release_version

        assert p.get_genome_download_link(
            "GRCh38.p13", mask=mask, toplevel=toplevel, version=version
        )

        sleep(1)

    @pytest.mark.skipif(travis and linux, reason="FTP does not work on Travis-Linux")
    def test_ensemblgenomes_genome_download_links(masking):
        """Test Ensembl FTP links for various genomes

        These genomes are hosted on ftp.ensemblgenomes.org.
        """
        p = genomepy.provider.ProviderBase.create("Ensembl")

        mask = masking if masking != "unmasked" else "none"

        for genome in ["Amel_4.5", "WBcel235"]:
            assert p.get_genome_download_link(genome, mask=mask)
            p.version = None  # reset version for next genome

        sleep(1)

    def test_ucsc_genome_download_links(masking):
        """Test UCSC HTTP links for various genomes

        Also test masking (unmasked should be ignored)."""
        p = genomepy.provider.ProviderBase.create("UCSC")

        for genome in ["sacCer3", "hg38"]:
            assert p.get_genome_download_link(genome, mask=masking)

        sleep(1)

    def test_ncbi_genome_download_links(masking):
        """Test NCBI HTTPS links for various genomes

        Also test masking (should be ignored).

        These genomes are hosted on ftp://ftp.ncbi.nlm.nih.gov."""
        p = genomepy.provider.ProviderBase.create("NCBI")

        for genome in ["Charlie1.0", "GRCh38.p13"]:
            assert p.get_genome_download_link(genome, mask=masking)

        sleep(1)
