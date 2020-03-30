import filecmp
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


def test_download_and_generate_annotation():
    out_dir = os.getcwd()
    localname = "my_annot"

    # unrecognized file type extension
    annot_url = "https://www.google.com"
    with pytest.raises(TypeError), TemporaryDirectory(dir=out_dir) as tmpdir:
        genomepy.provider.download_and_generate_annotation(
            genome_dir=tmpdir, annot_url=annot_url, localname=localname
        )


def test_attempt_download_and_report_back(capsys):
    out_dir = os.getcwd()
    localname = "my_annot"
    annot_url = "https://www.google.com"
    with pytest.raises(TypeError), TemporaryDirectory(dir=out_dir) as tmpdir:
        genomepy.provider.attempt_download_and_report_back(tmpdir, annot_url, localname)

    captured = capsys.readouterr().err.strip()
    assert captured.startswith(f"Using {annot_url}") and captured.endswith(
        "https://github.com/vanheeringen-lab/genomepy/issues"
    )


@pytest.mark.skipif(not travis, reason="Too slow")
def test_download_and_generate_annotation_and_attempt_download_and_report_back():
    out_dir = os.getcwd()
    localname = "my_annot"

    # fr3 from UCSC
    annot_url = "http://hgdownload.cse.ucsc.edu/goldenPath/fr3/database/ensGene.txt.gz"
    with TemporaryDirectory(dir=out_dir) as tmpdir:
        genomepy.provider.attempt_download_and_report_back(
            genome_dir=tmpdir, annot_url=annot_url, localname=localname
        )

        # check download_and_generate_annotation output
        fname = os.path.join(tmpdir, localname, localname + ".annotation.gtf.gz")
        validate_gzipped_gtf(fname)

        fname = os.path.join(tmpdir, localname, localname + ".annotation.bed.gz")
        validate_gzipped_bed(fname)

        # check attempt_download_and_report_back output
        readme = os.path.join(tmpdir, localname, "README.txt")
        with open(readme, "r") as f:
            assert f.readline() == f"Annotation url: {annot_url}\n"


def test_providerbase__init__():
    p = genomepy.provider.ProviderBase()

    result = sorted([x for x in dir(p) if not x.startswith("__")])
    expected = [
        "_providers",
        "create",
        "download_annotation",
        "download_genome",
        "get_genome_download_link",
        "list_install_options",
        "list_providers",
        "name",
        "register_provider",
        "tar_to_bigfile",
    ]

    assert result == expected
    assert p.name is None


def test_create():
    p = genomepy.provider.ProviderBase()
    p.create("Ensembl")

    with pytest.raises(ValueError):
        p.create("error")


def test_register_provider_and_list_providers():
    p = genomepy.provider.ProviderBase()
    assert isinstance(p._providers, dict)

    providers = ["ensembl", "ucsc", "ncbi", "url"]
    for provider in providers:
        assert provider in list(p.list_providers())


def test_list_install_options():
    p = genomepy.provider.ProviderBase()
    assert isinstance(p.list_install_options(), dict)
    assert len(p.list_install_options()) == 0

    with pytest.raises(ValueError):
        p.list_install_options(name="error")

    result = sorted(list(p.list_install_options(name="ensembl").keys()))
    expected = ["toplevel", "version"]
    assert result == expected


def test_tar_to_bigfile():
    p = genomepy.provider.ProviderBase()
    fname = "tests/data/tar.fa.tar.gz"
    outname = "tests/data/tar.fa"
    p.tar_to_bigfile(fname, outname)

    assert os.path.exists(outname)
    # tar.fa is a copy of gap.fa. Check if they are identical after untarring.
    assert filecmp.cmp(outname, "tests/data/gap.fa")
    os.unlink(outname)


def test_download_genome():
    pass
