import genomepy
import gzip
import os
import pytest
import shutil

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
    for field in p.taxid_fields + p.description_fields:
        assert field in genome
    assert genome["taxId"] == 9646


def test_assembly_accession(p):
    genome = p.genomes["sacCer3"]
    accession = p.assembly_accession(genome)

    assert accession.startswith("GCA_000146055")


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


# def test_list_available_genomes():
#     p = genomepy.provider.UcscProvider()
#     g = p.list_available_genomes()
#
#     assert "ailMel1" in g.keys()
#     assert g["ailMel1"]["scientificName"] == "Ailuropoda melanoleuca"
#     assert g["ailMel1"]["taxId"] == 9646
#
#
#
#
# def test_genome_taxid(name="sacCer3"):
#     p = genomepy.provider.UcscProvider()
#     t = p.genome_taxid(name)
#
#     assert t == 559292
#
#
# def test__genome_info_tuple():
#     p = genomepy.provider.UcscProvider()
#     t = p._genome_info_tuple("sacCer3")
#
#     assert isinstance(t, tuple)
#     assert t[0] == "sacCer3"
#     assert t[1] == "GCA_000146055.2"
#     assert t[2] == "Saccharomyces cerevisiae"
#     assert t[3] == "559292"
#
#
# def test_search():
#     p = genomepy.provider.UcscProvider()
#     for method in ["559292", "sacCer3"]:
#         s = p.search(method)
#         for genome in s:
#             if genome[0] == "sacCer3":
#                 assert genome[1] == "GCA_000146055.2"
#                 assert genome[3] == "559292"
#                 break
#
#
# def test_get_genome_download_link(name="sacCer3"):
#     p = genomepy.provider.UcscProvider()
#     link = p.get_genome_download_link(name)
#
#     assert link[0] == name
#     assert (
#         link[1]
#         == "http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz"
#     )
#
#
#
#
# @pytest.mark.skipif(not travis, reason="slow")
# def test_download_annotation(name="sacCer3"):
#     """Test UCSC annotation"""
#     p = genomepy.provider.UcscProvider()
#     out_dir = os.getcwd()
#     localname = "my_annot"
#
#     with TemporaryDirectory(dir=out_dir) as tmpdir:
#         p.download_annotation(name, tmpdir, localname=localname)
#
#         # check download_and_generate_annotation output
#         fname = os.path.join(tmpdir, localname, localname + ".annotation.gtf.gz")
#         validate_gzipped_gtf(fname)
#
#         fname = os.path.join(tmpdir, localname, localname + ".annotation.bed.gz")
#         validate_gzipped_bed(fname)
