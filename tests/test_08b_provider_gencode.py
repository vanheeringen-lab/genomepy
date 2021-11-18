import pytest

import genomepy
from tests import linux, travis


@pytest.mark.skipif(travis and linux, reason="FTP does not work on Travis-Linux")
def test_gencodeprovider(gencode):
    assert gencode.name == "GENCODE"
    assert gencode.taxid_fields == ["taxonomy_id"]


@pytest.mark.skipif(travis and linux, reason="FTP does not work on Travis-Linux")
def test_genome_info_tuple(gencode):
    t = gencode._genome_info_tuple("GRCh37")
    assert isinstance(t, tuple)
    assert t[0:4] == ("GRCh37", "GCA_000001405.1", 9606, True)


@pytest.mark.skipif(travis and linux, reason="FTP does not work on Travis-Linux")
def test_genomes(gencode):
    assert gencode.genomes["GRCh37"]["other_info"] == "GENCODE annotation + UCSC genome"
    assert gencode.genomes["GRCh38"]["assembly_accession"] == "GCA_000001405.15"


@pytest.mark.skipif(travis and linux, reason="FTP does not work on Travis-Linux")
def test_get_genome_download_link(gencode):
    link = gencode.get_genome_download_link("GRCh37", mask="soft")
    assert link in [
        "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz",
        "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.masked.gz",
    ]


@pytest.mark.skipif(travis and linux, reason="FTP does not work on Travis-Linux")
def test_get_annotation_download_links(gencode):
    # default annotation filing system
    genome = "GRCm39"
    annots = gencode.get_annotation_download_links(genome)
    expected = [  # release numbers removed
        "ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M",
        "/gencode.vM",
        ".annotation.gtf.gz",
    ]
    assert all([exp in annots[0] for exp in expected])

    # GRCh37, the one with the unique filing system.
    genome = "GRCh37"
    annots = gencode.get_annotation_download_links(genome)
    expected = [  # release numbers removed
        "ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_",
        "/GRCh37_mapping/gencode.v",
        "lift37.annotation.gtf.gz",
    ]
    assert all([exp in annots[0] for exp in expected])


@pytest.mark.skipif(travis and linux, reason="FTP does not work on Travis-Linux")
def test_download_annotation(gencode):
    gencode.download_annotation("GRCm39")  # smallest gencode annotation (0.8 GB)


def test_get_gencode2ucsc():
    genomes = {
        "test1": {"species": "Homo sapiens"},
        "test2": {"species": "Mus Musculus"},
        "test3": {"species": "whatever"},
    }
    gencode2ucsc = genomepy.providers.gencode.get_gencode2ucsc(genomes)
    assert gencode2ucsc["test1"] == "hg1"
    assert gencode2ucsc["test2"] == "mm2"
    assert gencode2ucsc["test3"] == "mm3"


def test_get_releases():
    listing = [
        "/path/to/release_22",
        "/path/to/mouse/release_M44",
        "/path/to/mouse/release_M33",
        "/path/to/release_01",  # too old
        "/path/to/something/else",
    ]
    specie = "human"
    releases = genomepy.providers.gencode.get_releases(listing, specie)

    assert releases == ["44", "33", "22"]

    specie = "mouse"
    releases = genomepy.providers.gencode.get_releases(listing, specie)

    assert releases == ["M44", "M33", "M22"]


def test_add_grch37():
    release = 42
    genomes = {
        "GRCh11": {},
        "GRCh22": {"annotations": [f"ftp/to/release_{release}/gtf"]},
    }
    genomes = genomepy.providers.gencode.add_grch37(genomes, "")
    expected = (
        f"/Gencode_human/release_{release}/GRCh37_mapping/"
        f"gencode.v{release}lift37.annotation.gtf.gz"
    )

    assert genomes["GRCh22"]["annotations"] == [f"ftp/to/release_{release}/gtf"]
    assert genomes["GRCh37"]["annotations"] == [expected]
