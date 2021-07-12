import os
from shutil import copyfile
from tempfile import TemporaryDirectory

import pytest

import genomepy.exceptions


def test_ucscprovider(ucsc):
    assert ucsc.name == "UCSC"
    assert ucsc.taxid_fields == ["taxId"]


def test__search_accession(ucsc):
    assert list(ucsc._search_accession("not_an_id")) == []
    assert next(ucsc._search_accession("GCA_000004335.1")) == "ailMel1"


def test_assembly_accession(ucsc):
    accession = ucsc.assembly_accession("sacCer3")

    assert accession.startswith("GCA_000146045")


def test_annotation_links(ucsc):
    annots = ucsc.annotation_links("sacCer3")
    expected = ["ensGene", "ncbiRefSeq"]
    assert annots == expected


def test_genome_info_tuple(ucsc):
    t = ucsc._genome_info_tuple("sacCer3")
    assert isinstance(t, tuple)
    assert t[0:4] == ("sacCer3", "GCA_000146045.2", 559292, [False, True, True, False])


def test_get_genome_download_link(ucsc):
    link = ucsc.get_genome_download_link("sacCer3", mask="soft")
    assert link in [
        "http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz",
        "http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz",
    ]

    link = ucsc.get_genome_download_link("danRer7", mask="hard")
    assert link in [
        "http://hgdownload.soe.ucsc.edu/goldenPath/danRer7/bigZips/chromFaMasked.tar.gz",
        "http://hgdownload.soe.ucsc.edu/goldenPath/danRer7/bigZips/danRer7.fa.masked.gz",
    ]


def test__post_process_download(ucsc):
    localname = "tmp"
    out_dir = os.getcwd()
    with TemporaryDirectory(dir=out_dir) as tmpdir:
        # this should skip without error
        ucsc._post_process_download(name=None, fname="", out_dir=None, mask="soft")
        ucsc._post_process_download(name=None, fname="", out_dir=None, mask="hard")
        ucsc._post_process_download(name=None, fname="", out_dir=None, mask="???")

        # copy fa file for unmasking
        g = os.path.join(tmpdir, localname + ".fa")
        copyfile("tests/data/gap.fa", g)

        ucsc._post_process_download(name=None, fname=g, out_dir=None, mask="none")
        assert os.path.exists(g)
        with open(g) as f:
            for line in f:
                assert "a" not in line


def test_get_annotation_download_links(ucsc):
    # any GTF format annotation
    genome = "sacCer3"
    annots = ucsc.get_annotation_download_links(genome)
    annots.sort()
    expected = ["ensGene", "ncbiRefSeq"]
    assert annots == expected

    # no annotation available
    assert ucsc.get_annotation_download_links("apiMel1") == []


def test_get_annotation_download_link(ucsc):
    # any annotation
    genome = "xenTro2"
    annot = ucsc.get_annotation_download_link(genome)
    assert annot == "ensGene"

    # specific annotation type
    genome = "sacCer3"
    annot = ucsc.get_annotation_download_link(
        genome, **{"ucsc_annotation_type": "NCBI_refseq"}
    )
    assert annot == "ncbiRefSeq"

    # no annotation available
    with pytest.raises(genomepy.exceptions.GenomeDownloadError):
        ucsc.get_annotation_download_link("apiMel1")

    # annotation type non-existing for this genome
    with pytest.raises(FileNotFoundError):
        ucsc.get_annotation_download_link(genome, **{"ucsc_annotation_type": "UCSC"})

    # annotation type non-existing
    with pytest.raises(FileNotFoundError):
        ucsc.get_annotation_download_link(genome, **{"ucsc_annotation_type": "what?"})


def test_get_genomes(ucsc):
    assert isinstance(ucsc.genomes, dict)
    assert "ailMel1" in ucsc.genomes
    genome = ucsc.genomes["ailMel1"]
    assert isinstance(genome, dict)
    for field in ucsc.accession_fields + ucsc.taxid_fields + ucsc.description_fields:
        assert field in genome
    assert genome["taxId"] == 9646
    assert genome["assembly_accession"] == "GCF_000004335.2"
    assert genome["annotations"] == ["ncbiRefSeq", "ensGene"]


def test_head_annotations(ucsc, caplog):
    ucsc.head_annotations("hg38", annotations=["ncbiRefSeq", "refGene"], n=1)

    assert "hg38_head_ncbiRefSeq.annotation.gtf" in caplog.text
    assert 'gene_name "DDX11L1";' in caplog.text

    assert "hg38_head_refGene.annotation.gtf" in caplog.text
    assert 'gene_name "WASH7P";' in caplog.text
