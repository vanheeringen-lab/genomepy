import genomepy
import gzip
import os
import pytest
import shutil

from tempfile import TemporaryDirectory


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
    return genomepy.provider.UcscProvider()


def test_ensemblprovider__init__(p):
    p2 = genomepy.provider.ProviderBase().create("UCSC")
    assert p.name == p2.name == "UCSC"
    assert p.taxid_fields == ["taxId"]


def test__get_genomes(p):
    assert isinstance(p.genomes, dict)
    assert "ailMel1" in p.genomes
    genome = p.genomes["ailMel1"]
    assert isinstance(genome, dict)
    for field in p.accession_fields + p.taxid_fields + p.description_fields:
        assert field in genome
    assert genome["taxId"] == 9646


def test_assembly_accession(p):
    genome = p.genomes["sacCer3"]
    accession = p.assembly_accession(genome)

    assert accession.startswith("GCA_000146045")


def test_genome_info_tuple(p):
    t = p._genome_info_tuple("sacCer3")
    assert isinstance(t, tuple)
    assert t[2:4] == ("Saccharomyces cerevisiae", "559292")


def test_get_genome_download_link(p):
    link = p.get_genome_download_link("sacCer3", mask="soft")
    assert (
        link
        == "http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz"
    )

    link = p.get_genome_download_link("danRer7", mask="hard")
    assert (
        link
        == "http://hgdownload.soe.ucsc.edu/goldenPath/danRer7/bigZips/danRer7.fa.masked.gz"
    )


def test__post_process_download(p):
    localname = "tmp"
    out_dir = os.getcwd()
    with TemporaryDirectory(dir=out_dir) as tmpdir:
        # this should skip without error
        p._post_process_download(
            name=None, localname=localname, out_dir=tmpdir, mask="soft"
        )
        p._post_process_download(
            name=None, localname=localname, out_dir=tmpdir, mask="hard"
        )
        p._post_process_download(
            name=None, localname=localname, out_dir=tmpdir, mask="???"
        )

        # copy fa file for unmasking
        g = os.path.join(tmpdir, localname + ".fa")
        shutil.copyfile("tests/data/gap.fa", g)

        p._post_process_download(
            name=None, localname=localname, out_dir=tmpdir, mask="none"
        )
        assert os.path.exists(g)
        with open(g) as f:
            for line in f:
                assert "a" not in line


def test_get_annotation_download_link(p):
    link = p.get_annotation_download_link("sacCer3")
    assert (
        link
        == "http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/database/ensGene.txt.gz"
    )
