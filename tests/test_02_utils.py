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


def test_get_genomes_dir(genomes_dir="tests/data"):
    # dir does exist
    gd = genomepy.utils.get_genomes_dir(genomes_dir=genomes_dir, check_exist=True)
    assert os.path.abspath(gd) == os.path.abspath(genomes_dir)

    # dir does not exist
    with pytest.raises(FileNotFoundError):
        genomepy.utils.get_genomes_dir(genomes_dir="fake_dir", check_exist=True)


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


def test_safe(unsafe_name="a name ", safe_name="a_name"):
    result = genomepy.utils.safe(unsafe_name)
    assert result == safe_name


def test_lower(unsafe_name="A nAme ", safe_name="a_name"):
    result = genomepy.utils.lower(unsafe_name)
    assert result == safe_name


def test_get_localname(name="XENTR_9.1", localname="my genome"):
    # localname input
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

    # Local path input
    path = "tests/data/sacCer3/sacCer3.fa"
    result = genomepy.utils.get_localname(name=path)
    assert result == "sacCer3"


def test_try_except_pass():
    def raise_error():
        raise ValueError

    # expected error is caught
    genomepy.utils.try_except_pass(ValueError, raise_error)
    genomepy.utils.try_except_pass((ValueError, FileNotFoundError), raise_error)

    # unexpected error is not caught
    with pytest.raises(ValueError):
        genomepy.utils.try_except_pass(FileNotFoundError, raise_error)
