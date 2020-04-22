import genomepy
import gzip
import os
import pytest

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
    return genomepy.provider.UrlProvider()


def test_ncbiprovider__init__(p):
    p2 = genomepy.provider.ProviderBase().create("URL")
    assert p2.name == p.name == "URL"
    assert p.genomes == {}


def test_genome_taxid(p):
    assert p.genome_taxid({}) == "na"


def test_assembly_accession(p):
    assert p.assembly_accession({}) == "na"


def test_search(p):
    result = p.search("???")
    expected = genomepy.provider.EnsemblProvider().search("???")
    assert isinstance(result, type(expected))
    with pytest.raises(StopIteration):
        assert next(result)


def test_list_install_options(p):
    result = list(p.list_install_options())
    expected = ["to_annotation"]
    assert result == expected


def test_get_genome_download_link(p):
    link = p.get_genome_download_link("url")
    assert link == "url"


def test_get_annotation_download_link(p):
    url = "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_Xenbase.gtf"
    link = p.get_annotation_download_link(None, **{"to_annotation": url})
    assert link == url

    with pytest.raises(TypeError):
        bad_url = "bad_url"
        p.get_annotation_download_link(None, **{"to_annotation": bad_url})

    with pytest.raises(TypeError):
        bad_url = "http://good_url.bad_ext"
        p.get_annotation_download_link(None, **{"to_annotation": bad_url})

    bad_url = "http://good_url_but_doesnt_exist.gtf"
    link = p.get_annotation_download_link(None, **{"to_annotation": bad_url})
    assert link is None


def test_search_url_for_annotation(p):
    url = "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XT9_1.fa.gz"
    link = p.search_url_for_annotation(url)
    expected = "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_Xenbase.gtf"
    assert link == expected


@pytest.mark.skipif(not travis or not linux, reason="slow")
def test_download_annotation(p):
    out_dir = os.getcwd()
    annot_url = "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_Xenbase.gtf"
    localname = "my_annot"
    with TemporaryDirectory(dir=out_dir) as tmpdir:
        p.download_annotation(
            url="string",
            genomes_dir=tmpdir,
            localname=localname,
            **{"to_annotation": annot_url},
        )

        # check download_and_generate_annotation output
        fname = os.path.join(tmpdir, localname, localname + ".annotation.gtf.gz")
        validate_gzipped_gtf(fname)

        fname = os.path.join(tmpdir, localname, localname + ".annotation.bed.gz")
        validate_gzipped_bed(fname)

        # check attempt_download_and_report_back output
        readme = os.path.join(tmpdir, localname, "README.txt")
        metadata, lines = genomepy.utils.read_readme(readme)
        assert metadata["annotation url"] == annot_url
