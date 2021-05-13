import os
from shutil import copyfile
from tempfile import TemporaryDirectory


def test_ucscprovider__init__(ucsc):
    assert ucsc.name == "UCSC"
    assert ucsc.taxid_fields == ["taxId"]


def test__get_genomes(ucsc):
    assert isinstance(ucsc.genomes, dict)
    assert "ailMel1" in ucsc.genomes
    genome = ucsc.genomes["ailMel1"]
    assert isinstance(genome, dict)
    for field in ucsc.accession_fields + ucsc.taxid_fields + ucsc.description_fields:
        assert field in genome
    assert genome["taxId"] == 9646


def test__search_accession(ucsc):
    assert list(ucsc._search_accession("not_an_id")) == []
    assert list(ucsc._search_accession("GCA_000004335.1")) == ["ailMel1"]


def test_assembly_accession(ucsc):
    genome = ucsc.genomes["sacCer3"]
    accession = ucsc.assembly_accession(genome)

    assert accession.startswith("GCA_000146045")


def test_genome_info_tuple(ucsc):
    t = ucsc._genome_info_tuple("sacCer3")
    assert isinstance(t, tuple)
    assert t[2:4] == ("Saccharomyces cerevisiae", "559292")


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


def test_get_annotation_download_link(ucsc):
    # any GTF format annotation
    genome = "sacCer3"
    link = ucsc.get_annotation_download_link(genome)
    assert genome in link
    assert link.endswith(".gtf.gz")

    # specific GTF annotation type
    link = ucsc.get_annotation_download_link(
        genome, **{"ucsc_annotation_type": "NCBI_refseq"}
    )
    assert genome in link
    assert link.endswith("ncbiRefSeq.gtf.gz")

    # any TXT format annotation (no GTF available)
    genome = "xenTro2"
    link = ucsc.get_annotation_download_link(genome)
    assert genome in link
    assert link.endswith(".txt.gz")

    # specific TXT annotation type (no GTF available)
    link = ucsc.get_annotation_download_link(
        genome, **{"ucsc_annotation_type": "UCSC_refseq"}
    )
    assert genome in link
    assert link.endswith("refGene.txt.gz")

    # non-existing annotation type
    link = ucsc.get_annotation_download_link(genome, **{"ucsc_annotation_type": "UCSC"})
    assert link is None
