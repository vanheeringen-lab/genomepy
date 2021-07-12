import filecmp
import os
from tempfile import TemporaryDirectory

import pytest

import genomepy.files
import genomepy.utils
from tests import travis
from tests.conftest import validate_annot


def test_baseprovider(base):
    assert base.name is None


def test__check_name(ucsc):
    ucsc._check_name("ailMel1")

    with pytest.raises(genomepy.exceptions.GenomeDownloadError):
        ucsc._check_name("not_a_real_genome")


def test__genome_info_tuple(base):
    with pytest.raises(NotImplementedError):
        base._genome_info_tuple(None)


def test_list_available_genomes(base):
    assert list(base.list_available_genomes()) == []


def test_genome_taxid(ucsc):
    taxid = ucsc.genome_taxid("ailMel1")
    assert taxid == 9646


def test_assembly_accession(ncbi):
    accession = ncbi.assembly_accession("AilMel_1.0")
    assert "000004335" in accession


def test_annotation_links(ncbi):
    # most genomes: 1 annotations
    links = ncbi.annotation_links("AilMel_1.0")
    expected = [
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/004/335/"
        "GCA_000004335.1_AilMel_1.0/GCA_000004335.1_AilMel_1.0_genomic.gff.gz"
    ]
    assert links == expected

    # some genomes: no annotation
    links = ncbi.annotation_links("Oryza_glaberrima_V1")
    assert links == []


def get_genome_download_link(base):
    with pytest.raises(NotImplementedError):
        base.get_genome_download_link()


@pytest.mark.skipif(not travis, reason="slow")
def test_download_genome(
    ucsc,
    name="sacCer3",
    localname="my_genome",
    mask="soft",
):
    out_dir = os.getcwd()

    with TemporaryDirectory(dir=out_dir) as tmpdir:
        ucsc.download_genome(
            name,
            genomes_dir=tmpdir,
            localname=localname,
            mask=mask,
        )

        genome = os.path.join(tmpdir, localname, localname + ".fa")
        assert os.path.exists(genome)

        readme = os.path.join(os.path.dirname(genome), "README.txt")
        with open(readme) as f:
            metadata = {}
            for line in f.readlines():
                vals = line.strip().split(":")
                metadata[vals[0].strip()] = (":".join(vals[1:])).strip()

        assert metadata["name"] == localname
        assert metadata["mask"] == mask


def test_get_annotation_download_links(base):
    with pytest.raises(NotImplementedError):
        base.get_annotation_download_links(None)


def test_get_annotation_download_link(ncbi):
    link = ncbi.get_annotation_download_link("AilMel_1.0")
    expected = (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/004/335/"
        "GCA_000004335.1_AilMel_1.0/GCA_000004335.1_AilMel_1.0_genomic.gff.gz"
    )
    assert link == expected


# @pytest.mark.skipif(not travis, reason="slow")
# def test__download_annotation(base):
#     out_dir = os.getcwd()
#     localname = "my_annot"
#
#     annot_url = "https://www.google.com"
#     with pytest.raises(TypeError), TemporaryDirectory(dir=out_dir) as tmpdir:
#         base._download_annotation(
#             genomes_dir=tmpdir, annot_url=annot_url, localname=localname
#         )
#
#     annot_url = "http://hgdownload.soe.ucsc.edu/goldenPath/ailMel1/bigZips/genes/ailMel1.ensGene.gtf.gz"
#     with TemporaryDirectory(dir=out_dir) as tmpdir:
#         base._download_annotation(
#             genomes_dir=tmpdir, annot_url=annot_url, localname=localname
#         )
#
#         fname = os.path.join(tmpdir, localname, localname + ".annotation.gtf")
#         validate_annot(fname, "gtf")
#
#         fname = os.path.join(tmpdir, localname, localname + ".annotation.bed")
#         validate_annot(fname, "bed")


@pytest.mark.skipif(not travis, reason="slow")
def test_download_annotation(ucsc):
    out_dir = os.getcwd()
    localname = "my_annot"

    name = "xenTro2"

    annot_urls = [
        "http://hgdownload.cse.ucsc.edu/goldenPath/xenTro2/database/ensGene.txt.gz",
        "http://hgdownload.cse.ucsc.edu/goldenPath/xenTro2/database/refGene.txt.gz",
    ]
    with TemporaryDirectory(dir=out_dir) as tmpdir:
        ucsc.download_annotation(name=name, genomes_dir=tmpdir, localname=localname)

        # check download_and_generate_annotation output
        fname = os.path.join(tmpdir, localname, localname + ".annotation.gtf")
        validate_annot(fname, "gtf")

        fname = os.path.join(tmpdir, localname, localname + ".annotation.bed")
        validate_annot(fname, "bed")

        # check attempt_download_and_report_back output
        readme = os.path.join(tmpdir, localname, "README.txt")
        metadata, lines = genomepy.files.read_readme(readme)
        assert metadata["annotation url"] in annot_urls


def test__search_text(ucsc):
    term = genomepy.utils.lower("Ailuropoda melanoleuca")
    assert list(ucsc._search_text("not_in_description")) == []
    assert next(ucsc._search_text(term)) == "ailMel1"


def test__search_accession(ncbi):
    assert list(ncbi._search_accession("not_an_id")) == []
    assert next(ncbi._search_accession("GCA_000004335.1")) == "AilMel_1.0"


def test__search_taxonomy(ucsc):
    assert list(ucsc._search_taxonomy("not_an_id")) == []
    assert next(ucsc._search_taxonomy("9646")) == "ailMel1"


def test_search(ucsc):
    for method in ["ailMel1", "9646", "Ailuropoda melanoleuca"]:
        genome = next(ucsc.search(method))
        assert genome[0] == "ailMel1"
        assert isinstance(genome, tuple)


def test_tar_to_bigfile():
    fname = "tests/data/tar2.fa.tar.gz"
    outname = "tests/data/tar2.fa"
    genomepy.providers.base.tar_to_bigfile(fname, outname)  # noqa

    assert os.path.exists(outname)
    # tar2.fa is a copy of tar1.fa. Check if they are identical after untarring.
    assert filecmp.cmp(outname, "tests/data/tar1.fa")
    os.unlink(outname)
