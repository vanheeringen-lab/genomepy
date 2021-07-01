import os
from shutil import copyfile
from tempfile import TemporaryDirectory

import pytest


def test_ucscprovider(ucsc):
    assert ucsc.name == "UCSC"
    assert ucsc.taxid_fields == ["taxId"]


def test__search_accession(ucsc):
    assert list(ucsc._search_accession("not_an_id")) == []
    assert next(ucsc._search_accession("GCA_000004335.1")) == "ailMel1"


def test_assembly_accession(ucsc):
    accession = ucsc.assembly_accession("sacCer3")

    assert accession.startswith("GCA_000146045")


def test__annot_types(ucsc):
    annots = ucsc._annot_types("sacCer3")
    expected = [False, True, True, False]
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
        ucsc._post_process_download(
            name=None, localname=localname, out_dir=tmpdir, mask="soft"
        )
        ucsc._post_process_download(
            name=None, localname=localname, out_dir=tmpdir, mask="hard"
        )
        ucsc._post_process_download(
            name=None, localname=localname, out_dir=tmpdir, mask="???"
        )

        # copy fa file for unmasking
        g = os.path.join(tmpdir, localname + ".fa")
        copyfile("tests/data/gap.fa", g)

        ucsc._post_process_download(
            name=None, localname=localname, out_dir=tmpdir, mask="none"
        )
        assert os.path.exists(g)
        with open(g) as f:
            for line in f:
                assert "a" not in line


def test_get_annotation_download_links(ucsc):
    # any GTF format annotation
    genome = "sacCer3"
    links = ucsc.get_annotation_download_links(genome)
    assert genome in links[0]
    assert links[0].endswith(".gtf.gz")

    # any TXT format annotation (no GTF available)
    genome = "xenTro2"
    links = ucsc.get_annotation_download_links(genome)
    assert genome in links[0]
    assert links[0].endswith(".txt.gz")

    # no annotation available
    assert ucsc.get_annotation_download_links("apiMel1") is None


def test_get_annotation_download_link(ucsc):
    # specific GTF annotation type
    genome = "sacCer3"
    link = ucsc.get_annotation_download_link(
        genome, **{"ucsc_annotation_type": "NCBI_refseq"}
    )
    assert genome in link
    assert link.endswith("ncbiRefSeq.gtf.gz")

    # specific TXT annotation type (no GTF available)
    genome = "xenTro2"
    link = ucsc.get_annotation_download_link(
        genome, **{"ucsc_annotation_type": "UCSC_refseq"}
    )
    assert genome in link
    assert link.endswith("refGene.txt.gz")

    # annotation type non-existing for this genome
    with pytest.raises(FileNotFoundError):
        ucsc.get_annotation_download_link(genome, **{"ucsc_annotation_type": "UCSC"})

    # annotation type non-existing
    with pytest.raises(ValueError):
        ucsc.get_annotation_download_link(genome, **{"ucsc_annotation_type": "what?"})

    # no annotation available
    with pytest.raises(FileNotFoundError):
        ucsc.get_annotation_download_link("apiMel1")


def test_get_genomes(ucsc):
    assert isinstance(ucsc.genomes, dict)
    assert "ailMel1" in ucsc.genomes
    genome = ucsc.genomes["ailMel1"]
    assert isinstance(genome, dict)
    for field in ucsc.accession_fields + ucsc.taxid_fields + ucsc.description_fields:
        assert field in genome
    assert genome["taxId"] == 9646
