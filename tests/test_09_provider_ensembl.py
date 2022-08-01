import pytest
import requests

import genomepy.providers.ensembl


def test_ensemblprovider(ensembl):
    assert ensembl.name == "Ensembl"
    assert ensembl.taxid_fields == ["taxonomy_id"]


def test_genome_info_tuple(ensembl):
    t = ensembl._genome_info_tuple("KH", size=True)
    assert isinstance(t, tuple)
    assert t[:-1] == (
        "KH",
        "GCA_000224145.1",
        7719,
        True,
        "Ciona intestinalis",
        115227500,
    )


def test_get_version(ensembl):
    v = ensembl.get_version(vertebrates=True)
    assert v.isnumeric()
    assert int(v) > 100

    v = ensembl.get_version()
    assert v.isnumeric()
    assert int(v) > 50

    v = ensembl.get_version(set_version=42)
    assert v.isnumeric()
    assert v == "42"


def test_get_division(ensembl):
    genome = ensembl.genomes["TAIR10"]
    division, is_vertebrate = genomepy.providers.ensembl.get_division(genome)
    assert division == "plants"
    assert is_vertebrate is False

    genome = ensembl.genomes["GRCh38.p13"]
    division, is_vertebrate = genomepy.providers.ensembl.get_division(genome)
    assert division == "vertebrates"
    assert is_vertebrate is True


def test_get_genome_download_link(ensembl):
    # non vertebrate: soft masked
    link = ensembl.get_genome_download_link("TAIR10", mask="soft", **{"version": 46})
    assert (
        link
        == "http://ftp.ensemblgenomes.org/pub/release-46/plants/"
        + "fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz"
    )

    # non vertebrate: < release 30
    link = ensembl.get_genome_download_link("TAIR10", mask="soft", **{"version": 25})
    assert (
        link
        == "http://ftp.ensemblgenomes.org/pub/release-25/plants/"
        + "fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.25.dna_sm.toplevel.fa.gz"
    )

    # non vertebrate: entry not using url_name (see issue #205)
    link = ensembl.get_genome_download_link("WBcel235", mask="soft", **{"version": 35})
    assert (
        link
        == "http://ftp.ensemblgenomes.org/pub/release-35/metazoa/"
        + "fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa.gz"
    )

    # vertebrate with primary assembly: unmasked
    link = ensembl.get_genome_download_link("GRCz11", mask="none", **{"version": 98})
    assert (
        link
        == "http://ftp.ensembl.org/pub/release-98/fasta/"
        + "danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz"
    )

    # vertebrate with primary assembly: hard masked and toplevel only
    link = ensembl.get_genome_download_link(
        "GRCz11", mask="hard", **{"version": 98, "toplevel": True}
    )
    assert (
        link
        == "http://ftp.ensembl.org/pub/release-98/fasta/"
        + "danio_rerio/dna/Danio_rerio.GRCz11.dna_rm.toplevel.fa.gz"
    )

    # vertebrate: latest version
    version = ensembl.get_version(True)
    link = ensembl.get_genome_download_link("GRCz11", **{"toplevel": True})
    expected_link = (
        f"http://ftp.ensembl.org/pub/release-{version}/"
        "fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna_sm.toplevel.fa.gz"
    )
    assert link == expected_link


def test_get_annotation_download_links(ensembl):
    # non vertebrate
    links = ensembl.get_annotation_download_links("TAIR10", **{"version": 46})
    expected_link = (
        "http://ftp.ensemblgenomes.org/pub/release-46/plants/"
        "gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.46.gtf.gz"
    )
    assert links[0] == expected_link

    # vertebrate
    links = ensembl.get_annotation_download_links("GRCz11", **{"version": 98})
    expected_link = (
        "http://ftp.ensembl.org/pub/release-98/gtf/"
        "danio_rerio/Danio_rerio.GRCz11.98.gtf.gz"
    )
    assert links[0] == expected_link

    # non vertebrate: latest version
    version = ensembl.get_version(False)
    links = ensembl.get_annotation_download_links("TAIR10")
    expected_link = (
        f"http://ftp.ensemblgenomes.org/pub/release-{version}/plants/"
        f"gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.{version}.gtf.gz"
    )
    assert links[0] == expected_link

    # vertebrate: latest version
    version = ensembl.get_version(True)
    links = ensembl.get_annotation_download_links("GRCz11")
    expected_link = (
        f"http://ftp.ensembl.org/pub/release-{version}/"
        f"gtf/danio_rerio/Danio_rerio.GRCz11.{version}.gtf.gz"
    )
    assert links[0] == expected_link


def test_request_json():
    divisions = genomepy.providers.ensembl.request_json(
        "https://rest.ensembl.org/", "info/divisions?"
    )
    assert isinstance(divisions, list)
    assert "EnsemblVertebrates" in divisions

    # test not r.ok
    with pytest.raises(requests.exceptions.HTTPError):
        genomepy.providers.ensembl.request_json("https://rest.ensembl.org/", "error")


def test_get_genomes(ensembl):
    assert isinstance(ensembl.genomes, dict)
    assert "KH" in ensembl.genomes
    genome = ensembl.genomes["KH"]
    assert isinstance(genome, dict)
    for field in (
        ensembl.accession_fields + ensembl.taxid_fields + ensembl.description_fields
    ):
        assert field in genome
    assert genome["taxonomy_id"] == 7719
