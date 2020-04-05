import genomepy
import gzip
import os
import pytest
import requests

from tempfile import TemporaryDirectory
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


def test_ensemblprovider__init__():
    p = genomepy.provider.EnsemblProvider()
    p2 = genomepy.provider.ProviderBase().create("Ensembl")
    assert p2.name == p.name == "Ensembl"


def test_request_json():
    p = genomepy.provider.EnsemblProvider()
    divisions = p.request_json("info/divisions?")
    assert isinstance(divisions, list)
    assert "EnsemblVertebrates" in divisions

    # test not r.ok
    with pytest.raises(requests.exceptions.HTTPError):
        p.request_json("error")


def test_ensembl_list_install_options():
    p = genomepy.provider.EnsemblProvider()
    result = sorted(list(p.list_install_options(name="ensembl").keys()))
    expected = ["toplevel", "version"]
    assert result == expected


def test_ensembl_list_available_genomes():
    p = genomepy.provider.EnsemblProvider()

    # test dict for essential terms
    g = p.list_available_genomes(as_dict=True)
    g1 = next(g)
    result = sorted(list(g1.keys()))
    expected = [
        "assembly_accession",
        "assembly_name",
        "division",
        "genebuild",
        "name",
        "scientific_name",
        "taxonomy_id",
        "url_name",
    ]

    assert isinstance(g1, dict)
    for key in expected:
        assert key in result

    # test tuple
    g = p.list_available_genomes(as_dict=False)
    for genome in g:
        assert isinstance(genome, tuple)
        # for the (currently) first genome in the list: check for matching output
        if genome[0] == "ENA_1":
            assert genome[1] == "albugo_laibachii"
            break


def test_ensembl__get_genome_info(name="KH"):
    p = genomepy.provider.EnsemblProvider()
    g = p._get_genome_info(name)

    assert g["assembly_name"] == name
    assert g["assembly_accession"] == "GCA_000224145.1"
    assert g["taxonomy_id"] == 7719

    with pytest.raises(genomepy.exceptions.GenomeDownloadError):
        p._get_genome_info("error")


def test_ensembl_assembly_accession(name="KH"):
    p = genomepy.provider.EnsemblProvider()
    a = p.assembly_accession(name)

    assert a == "GCA_000224145.1"


def test_ensembl_genome_taxid(name="KH"):
    p = genomepy.provider.EnsemblProvider()
    t = p.genome_taxid(name)

    assert t == 7719


def test_ensembl__genome_info_tuple():
    p = genomepy.provider.EnsemblProvider()
    g = p.list_available_genomes(as_dict=True)

    for genome in g:
        t = p._genome_info_tuple(genome)
        assert isinstance(t, tuple)
        # for the (currently) first genome in the list: check for matching output
        if t[0] == "ENA_1":
            assert t[2] == "Albugo laibachii"
            break


def test_ensembl_search():
    p = genomepy.provider.EnsemblProvider()
    for method in ["7719", "KH"]:
        s = p.search(method)
        for genome in s:
            if genome[0] == "KH":
                assert genome[1] == "GCA_000224145.1"
                assert genome[3] == "7719"
                break


def test_ensembl_get_version():
    p = genomepy.provider.EnsemblProvider()
    ftp_site = "http://ftp.ensembl.org/pub"
    v = p.get_version(ftp_site)

    # note: this test will break every time ensembl releases a new version
    assert v == "99"


def test_ensembl_get_genome_download_link(name="KH"):
    p = genomepy.provider.EnsemblProvider()
    link = p.get_genome_download_link(name)

    assert link[0] == name
    assert (
        link[1]
        == "http://ftp.ensembl.org/pub//release-99/fasta/ciona_intestinalis"
        + "/dna//Ciona_intestinalis.KH.dna_sm.toplevel.fa.gz"
    )


@pytest.mark.skipif(not travis, reason="slow")
def test_ensembl_download_annotation(name="KH", version=98):
    """Test Ensembl annotation

    This annotation is hosted on https://ftp.ensembl.org.
    """
    p = genomepy.provider.EnsemblProvider()
    out_dir = os.getcwd()
    localname = "my_annot"

    # Only test on vertebrates as these are downloaded over HTTP.
    # All others are downloaded over FTP, which is unreliable on Travis.
    with TemporaryDirectory(dir=out_dir) as tmpdir:
        p.download_annotation(name, tmpdir, localname=localname, version=version)

        # check download_and_generate_annotation output
        fname = os.path.join(tmpdir, localname, localname + ".annotation.gtf.gz")
        validate_gzipped_gtf(fname)

        fname = os.path.join(tmpdir, localname, localname + ".annotation.bed.gz")
        validate_gzipped_bed(fname)
