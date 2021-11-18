import pytest


def test_localprovider(local):
    assert local.name == "Local"
    assert local.genomes == {}


def test_genome_taxid(local):
    assert local.genome_taxid(None) is None


def test_assembly_accession(local):
    assert local.assembly_accession(None) is None


def test_search(local):
    result = local.search("???")
    with pytest.raises(StopIteration):
        assert next(result)


def test__genome_info_tuple(local):
    assert local._genome_info_tuple(None) is tuple()


def test__check_name(local):
    assert local._check_name(None) is None


def test_get_genome_download_link(local):
    target = "tests/data/sacCer3/sacCer3.fa"
    link = local.get_genome_download_link(target)
    assert link.endswith(target)

    with pytest.raises(FileNotFoundError):
        local.get_genome_download_link("this path does not exist")


def test_get_annotation_download_link(local):
    target = "tests/data/sacCer3/sacCer3.annotation.gtf"
    link = local.get_annotation_download_link(None, **{"path_to_annotation": target})
    assert link.endswith(target)

    with pytest.raises(FileNotFoundError):
        bad_path = "bad/path"
        local.get_annotation_download_link(None, **{"path_to_annotation": bad_path})

    with pytest.raises(TypeError):
        bad_ext = "tests/data/sacCer3/sacCer3.fa"
        local.get_annotation_download_link(None, **{"path_to_annotation": bad_ext})


def test_get_annotation_download_links(local):
    target = "tests/data/sacCer3/sacCer3.fa"
    expected = "tests/data/sacCer3/sacCer3.annotation.gtf"
    links = local.get_annotation_download_links(target)
    assert links[0].endswith(expected)
