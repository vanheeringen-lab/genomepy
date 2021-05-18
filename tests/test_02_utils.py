import filecmp
import os
from tempfile import NamedTemporaryFile
import shutil
import subprocess as sp
import urllib.request
import urllib.error

import pytest

import genomepy
import genomepy.utils
from tests import linux


def test_read_readme():
    metadata, lines = genomepy.utils.read_readme("tests/data/not_an_existsing_file")
    assert metadata["name"] == "na"
    assert metadata["sanitized annotation"] == "no"
    assert lines == []

    wd = os.getcwd()
    readme = os.path.join(wd, "README.txt")
    with open(readme, "w") as f:
        f.writelines("provider: asd\n")
        f.writelines(
            "this is a regular line\n  \n  \n \nonly one empty line before me\n"
        )

    metadata, lines = genomepy.utils.read_readme(readme)
    assert metadata["provider"] == "asd"
    assert lines == ["this is a regular line", "", "only one empty line before me"]

    os.unlink(readme)


def test_write_readme():
    wd = os.getcwd()
    readme = os.path.join(wd, "README.txt")
    metadata, lines = genomepy.utils.read_readme(readme)

    metadata["name"] = "my_cool_genome"
    lines = ["", "I wanted to do some regex, but regex is hard"]
    genomepy.utils.write_readme(readme, metadata, lines)

    metadata2, lines2 = genomepy.utils.read_readme(readme)
    assert metadata == metadata2
    assert lines == lines2

    os.unlink(readme)


def test_update_readme():
    wd = os.getcwd()
    readme = os.path.join(wd, "README.txt")

    # new data added
    updated_metadata = {"key": "value"}
    extra_lines = ["let's add some text to the lines section"]
    genomepy.utils.update_readme(readme, updated_metadata, extra_lines)
    metadata, lines = genomepy.utils.read_readme(readme)
    assert metadata["key"] == updated_metadata["key"]
    assert extra_lines[0] in lines

    # new data overwrites
    updated_metadata = {"key": "value2"}
    genomepy.utils.update_readme(readme, updated_metadata)
    metadata, lines = genomepy.utils.read_readme(readme)
    assert metadata["key"] == updated_metadata["key"]

    os.unlink(readme)


def test_generate_gap_bed(fname="tests/data/gap.fa", outname="tests/data/gap.bed"):
    tmp = NamedTemporaryFile().name
    genomepy.utils.generate_gap_bed(fname, tmp)

    result = open(tmp).read()
    expected = open(outname).read()

    assert result == expected


def test_generate_fa_sizes(infa="tests/data/gap.fa"):
    tmp = NamedTemporaryFile().name
    genomepy.utils.generate_fa_sizes(infa, tmp)

    result = open(tmp).read()

    assert result == "chr1\t28\nchr2\t45\nchr3\t15\n"


def test_filter_fasta(fname="tests/data/regexp/regexp.fa"):
    # no sequences left after filtering
    with pytest.raises(Exception), NamedTemporaryFile(suffix=".fa") as tmpfa:
        regex = "missing_chromosome"
        genomepy.utils.filter_fasta(fname, regex=regex, outfa=tmpfa.name)

    # function proper
    regexps = [
        ("Chr.*", 2, 15),
        ("Scaffold.*", 1, 16),
        ("scaffold_.*", 3, 14),
        (r"^\d+$", 4, 13),
        ("chr.*", 4, 13),
    ]
    for regex, match, no_match in regexps:
        with NamedTemporaryFile(suffix=".fa") as tmpfa:
            keys = genomepy.utils.filter_fasta(
                fname, regex=regex, invert_match=False, outfa=tmpfa.name
            ).keys()
            assert len(keys) == match, regex

        with NamedTemporaryFile(suffix=".fa") as tmpfa:
            keys = genomepy.utils.filter_fasta(
                fname, regex=regex, invert_match=True, outfa=tmpfa.name
            ).keys()
            assert len(keys) == no_match, regex


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


def test_glob_ext_files(file="tests/data/small_genome.fa"):
    assert file not in genomepy.utils.glob_ext_files("tests/data")
    assert file + ".gz" in genomepy.utils.glob_ext_files("tests/data")
    assert len(genomepy.utils.glob_ext_files("tests/data", "fake_ext")) == 0
    assert len(genomepy.utils.glob_ext_files("tests/data/regexp")) == 1


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


def test_tar_to_bigfile():
    fname = "tests/data/tar2.fa.tar.gz"
    outname = "tests/data/tar2.fa"
    genomepy.utils.tar_to_bigfile(fname, outname)

    assert os.path.exists(outname)
    # tar2.fa is a copy of tar1.fa. Check if they are identical after untarring.
    assert filecmp.cmp(outname, "tests/data/tar1.fa")
    os.unlink(outname)


def test_gunzip_and_name(fname="tests/data/small_genome.fa.gz"):
    assert os.path.exists(fname)
    fname, gzip_file = genomepy.utils.gunzip_and_name(fname)
    assert gzip_file and fname.endswith(".fa")
    assert os.path.exists(fname)
    assert not os.path.exists(fname + ".gz")


def test_gzip_and_name(fname="tests/data/small_genome.fa"):
    assert os.path.exists(fname)
    fname = genomepy.utils.gzip_and_name(fname)
    assert fname.endswith(".gz")
    assert os.path.exists(fname)
    assert not os.path.exists(fname[:-3])

    fname, _ = genomepy.utils.gunzip_and_name(fname)
    assert fname.endswith(".fa")
    assert os.path.exists(fname)
    assert not os.path.exists(fname + ".gz")


def test_bgzip_and_name(fname="tests/data/small_genome.fa"):
    assert os.path.exists(fname)
    fname = genomepy.utils.bgzip_and_name(fname)
    assert fname.endswith(".gz")
    assert os.path.exists(fname)
    assert not os.path.exists(fname[:-3])

    with pytest.raises(sp.CalledProcessError):
        genomepy.utils.bgzip_and_name("tests/data/nofile.fa")


def test_retry(capsys):
    # runs passed function
    txt = "hello world"
    genomepy.utils.retry(print, 1, txt)
    captured = capsys.readouterr().out.strip()
    assert captured == txt

    # handles URLErrors
    def _offline_func():
        raise urllib.error.URLError("this function is offline")

    assert genomepy.utils.retry(_offline_func, 1) is None


def test_check_url():
    assert genomepy.utils.check_url("http://ftp.xenbase.org/pub/Genomics/JGI/README")
    # wrong/offline urls:
    assert not genomepy.utils.check_url("https://www.thiswebsiteisoffline.nl/")


def test_read_url(
    url="http://ftp.xenbase.org/pub/Genomics/JGI/README", expected="The data"
):
    text = genomepy.utils.read_url(url)
    assert text.startswith(expected)


def test_get_file_info(fname="tests/data/small_genome.fa.gz"):
    ext, gz = genomepy.utils.get_file_info(fname)
    assert ext == ".fa" and gz

    ext, gz = genomepy.utils.get_file_info(fname[:-2] + "fai")
    assert ext == ".fai" and not gz


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
