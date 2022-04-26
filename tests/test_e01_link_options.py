import pytest

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

    def test_ensembl_genome_download_links(assembly, masking, release_version, ensembl):
        """Test Ensembl links with various options

        These genomes are hosted on ftp.ensembl.org

        Vertebrates are downloaded from HTTP.
        """
        mask = masking if masking != "unmasked" else "none"
        toplevel = False if assembly == "primary_assembly" else True
        version = release_version
        assert ensembl.get_genome_download_link(
            "GRCh38.p13", mask=mask, toplevel=toplevel, version=version
        )

    def test_ensemblgenomes_genome_download_links(masking, ensembl):
        """Test Ensembl FTP links for various genomes

        These genomes are hosted on ftp.ensemblgenomes.org.
        """
        mask = masking if masking != "unmasked" else "none"
        for genome in ["Amel_HAv3.1", "ASM23943v1"]:
            assert ensembl.get_genome_download_link(genome, mask=mask)

    def test_ucsc_genome_download_links(masking, ucsc):
        """Test UCSC HTTP links for various genomes

        Also test masking (unmasked should be ignored)."""
        for genome in ["sacCer3", "hg38"]:
            assert ucsc.get_genome_download_link(genome, mask=masking)

    def test_ncbi_genome_download_links(masking, ncbi):
        """Test NCBI HTTPS links for various genomes

        Also test masking (should be ignored).

        These genomes are hosted on ftp://ftp.ncbi.nlm.nih.gov."""
        for genome in ["Charlie1.0", "GRCh38.p13"]:
            assert ncbi.get_genome_download_link(genome, mask=masking)
