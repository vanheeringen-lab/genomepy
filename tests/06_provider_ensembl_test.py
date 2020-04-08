import genomepy
import gzip
import os
import pytest
import requests

from platform import system

linux = system() == "Linux"
travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"


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


@pytest.fixture(scope="module")
def p():
    return genomepy.provider.EnsemblProvider()


def test_ensemblprovider__init__(p):
    p2 = genomepy.provider.ProviderBase().create("Ensembl")
    assert p.name == p2.name == "Ensembl"
    assert p.taxid_fields == ["taxonomy_id"]


def test__request_json(p):
    divisions = p._request_json("info/divisions?")
    assert isinstance(divisions, list)
    assert "EnsemblVertebrates" in divisions

    # test not r.ok
    with pytest.raises(requests.exceptions.HTTPError):
        p._request_json("error")


def test__get_genomes(p):
    assert isinstance(p.genomes, dict)
    assert "KH" in p.genomes
    genome = p.genomes["KH"]
    assert isinstance(genome, dict)
    for field in p.accession_fields + p.taxid_fields + p.description_fields:
        assert field in genome
    assert genome["taxonomy_id"] == 7719


def test_list_install_options(p):
    result = sorted(list(p.list_install_options(name="ensembl").keys()))
    expected = ["toplevel", "version"]
    assert result == expected


def test_genome_info_tuple(p):
    t = p._genome_info_tuple("KH")
    assert isinstance(t, tuple)
    assert t[2:4] == ("Ciona intestinalis", "7719")


def test_get_version(p):
    # note: this test will break every time ensembl releases a new version
    ftp_site = "http://ftp.ensembl.org/pub"
    v = p.get_version(ftp_site)
    assert v == "99"

    ftp_site = "ftp://ftp.ensemblgenomes.org/pub"
    v = p.get_version(ftp_site)
    assert v == "46"


def test_get_genome_download_link(p):
    # non vertebrate: soft masked
    link = p.get_genome_download_link("TAIR10", mask="soft", **{"version": 46})
    assert (
        link
        == "ftp://ftp.ensemblgenomes.org/pub/plants/release-46/"
        + "fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz"
    )

    # vertebrate with primary assembly: unmasked
    link = p.get_genome_download_link("GRCz11", mask="none", **{"version": 98})
    assert (
        link
        == "http://ftp.ensembl.org/pub/release-98/fasta/"
        + "danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz"
    )

    # vertebrate with primary assembly: hard masked and toplevel only
    link = p.get_genome_download_link(
        "GRCz11", mask="hard", **{"version": 98, "toplevel": True}
    )
    assert (
        link
        == "http://ftp.ensembl.org/pub/release-98/fasta/"
        + "danio_rerio/dna/Danio_rerio.GRCz11.dna_rm.toplevel.fa.gz"
    )

    # vertebrate: any version
    p.version = None
    link = p.get_genome_download_link("GRCz11", **{"toplevel": True})
    for substring in [
        "http://ftp.ensembl.org/pub/release-",
        "/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna_sm.toplevel.fa.gz",
    ]:
        assert substring in link


def test_get_annotation_download_link(p):
    # non vertebrate
    link = p.get_annotation_download_link("TAIR10", **{"version": 46})
    assert (
        link
        == "ftp://ftp.ensemblgenomes.org/pub/plants/release-46/"
        + "gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.46.gtf.gz"
    )

    # vertebrate
    link = p.get_annotation_download_link("GRCz11", **{"version": 98})
    assert (
        link
        == "http://ftp.ensembl.org/pub/release-98/gtf/"
        + "danio_rerio/Danio_rerio.GRCz11.98.gtf.gz"
    )

    # vertebrate: any version
    p.version = None
    link = p.get_annotation_download_link("GRCz11")
    for substring in [
        "http://ftp.ensembl.org/pub/release-",
        "/gtf/danio_rerio/Danio_rerio.GRCz11.",
        ".gtf.gz",
    ]:
        assert substring in link
