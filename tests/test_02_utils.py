import os
import shutil

import pytest

import genomepy.utils
from tests import linux


def test_mkdir_p(path="./tests/dir1/dir2/nestled_dir"):
    shutil.rmtree("./tests/dir1", ignore_errors=True)
    assert not os.path.isdir(path)

    # make a non-existing nestled-dir
    genomepy.utils.mkdir_p(path)
    assert os.path.isdir(path)

    # try to make an existing dir
    genomepy.utils.mkdir_p(path)

    shutil.rmtree("./tests/dir1", ignore_errors=True)
    assert not os.path.isdir(path)


def test_rm_rf(path="./tests/dir1/dir2/nestled_dir"):
    genomepy.utils.mkdir_p(path)
    fname = os.path.join(path, "file.txt")
    with open(fname, "w") as f:
        f.write("hello world!")

    # try to remove an existing file
    assert os.path.isfile(fname)
    genomepy.utils.rm_rf(fname)
    assert not os.path.isfile(fname)

    # try to remove an existing dir
    assert os.path.isdir(path)
    genomepy.utils.rm_rf(path)
    assert not os.path.isdir(path)

    # try to remove a non-existing dir
    genomepy.utils.rm_rf(path)


def test_cmd_ok():
    assert genomepy.utils.cmd_ok("STAR")
    assert genomepy.utils.cmd_ok("bwa")
    assert not genomepy.utils.cmd_ok("missing_cmd")


def test_run_index_cmd(caplog, name="tests", good_cmd="ls", bad_cmd="bad_cmd"):
    genomepy.utils.run_index_cmd(name=name, cmd=good_cmd)

    # bad_command not found error
    genomepy.utils.run_index_cmd(name=name, cmd=bad_cmd)
    # captured = capsys.readouterr().err.strip()
    if linux:
        assert f"{bad_cmd}: not found" in caplog.text
    else:
        # thanks for being different mac. That took me 30 minutes of my life...
        assert f"{bad_cmd}: command not found" in caplog.text


def test_get_genomes_dir(genomes_dir="tests/data"):
    # dir does exist
    gd = genomepy.utils.get_genomes_dir(genomes_dir=genomes_dir, check_exist=True)
    assert os.path.abspath(gd) == os.path.abspath(genomes_dir)

    # dir does not exist
    with pytest.raises(FileNotFoundError):
        genomepy.utils.get_genomes_dir(genomes_dir="fake_dir", check_exist=True)


def test_safe(unsafe_name="a name ", safe_name="a_name"):
    result = genomepy.utils.safe(unsafe_name)
    assert result == safe_name


def test_get_localname(name="XENTR_9.1", localname="my genome"):
    # name + localname input
    result = genomepy.utils.get_localname(name=name, localname=localname)
    assert result == genomepy.utils.safe(localname)

    # name input
    result = genomepy.utils.get_localname(name=name)
    assert result == genomepy.utils.safe(name)

    # URL input (simple)
    url = "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_genome.fa.gz"
    result = genomepy.utils.get_localname(name=url)
    assert result == name

    # URL input (complex)
    url2 = (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/"
        "GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz"
    )
    result = genomepy.utils.get_localname(name=url2)
    assert result == "GCF_000027325.1_ASM2732v1"


def test__open():
    # read/write regular file
    a1 = "tests/data/data.annotation.gtf"
    with genomepy.annotation._open(a1, "w") as gtf:
        gtf.write("regular file")

    with open(a1) as gtf:
        lines1 = gtf.readlines()
    assert lines1 == ["regular file"]

    with genomepy.annotation._open(a1) as gtf:
        lines2 = gtf.readlines()
    assert lines2 == lines1
    genomepy.utils.rm_rf(a1)

    # read/write gzipped file
    a2 = "tests/data/data.annotation.gtf.gz"
    with genomepy.annotation._open(a2, "w") as gtf:
        gtf.write("gzipped file")

    with pytest.raises(UnicodeDecodeError):
        with open(a2) as gtf:
            gtf.read()

    with genomepy.annotation._open(a2) as gtf:
        lines = gtf.readlines()
    assert lines == ["gzipped file"]
    genomepy.utils.rm_rf(a1)


def test_best_search_result():
    # result of list(genomepy.ProviderBase.search_all("GCA_000001405", provider="NCBI"))
    ncbi_human = [
        [
            "GRCh38",
            "NCBI",
            "GCF_000001405.26",
            "Homo sapiens",
            "9606",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p1",
            "NCBI",
            "GCF_000001405.27",
            "Homo sapiens",
            "9606",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p2",
            "NCBI",
            "GCF_000001405.28",
            "Homo sapiens",
            "9606",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p3",
            "NCBI",
            "GCF_000001405.29",
            "Homo sapiens",
            "9606",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p4",
            "NCBI",
            "GCF_000001405.30",
            "Homo sapiens",
            "9606",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p5",
            "NCBI",
            "GCF_000001405.31",
            "Homo sapiens",
            "9606",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p6",
            "NCBI",
            "GCF_000001405.32",
            "Homo sapiens",
            "9606",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p7",
            "NCBI",
            "GCF_000001405.33",
            "Homo sapiens",
            "9606",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p8",
            "NCBI",
            "GCF_000001405.34",
            "Homo sapiens",
            "9606",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p9",
            "NCBI",
            "GCF_000001405.35",
            "Homo sapiens",
            "9606",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p10",
            "NCBI",
            "GCF_000001405.36",
            "Homo sapiens",
            "9606",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p11",
            "NCBI",
            "GCF_000001405.37",
            "Homo sapiens",
            "9606",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p12",
            "NCBI",
            "GCF_000001405.38",
            "Homo sapiens",
            "9606",
            "Genome Reference Consortium",
        ],
        [
            "GRCh37.p13",
            "NCBI",
            "GCF_000001405.25",
            "Homo sapiens",
            "9606",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p13",
            "NCBI",
            "GCF_000001405.39",
            "Homo sapiens",
            "9606",
            "Genome Reference Consortium",
        ],
    ]

    # list of search results
    best_result = genomepy.utils.best_search_result("GCA_000001405.39", ncbi_human)
    assert best_result[2] == "GCF_000001405.39"

    # zero search results
    assert len(genomepy.utils.best_search_result("GCA_000001405.39", [])) == 0

    # one search results
    best_result = genomepy.utils.best_search_result("GCA_000001405.26", [ncbi_human[0]])
    assert best_result[2] == ncbi_human[0][2]
