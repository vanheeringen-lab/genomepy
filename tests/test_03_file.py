import filecmp
import os
import subprocess as sp
from tempfile import NamedTemporaryFile

import pytest

import genomepy.files


def test_read_readme():
    metadata, lines = genomepy.files.read_readme("tests/data/not_an_existsing_file")
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

    metadata, lines = genomepy.files.read_readme(readme)
    assert metadata["provider"] == "asd"
    assert lines == ["this is a regular line", "", "only one empty line before me"]

    os.unlink(readme)


def test_write_readme():
    wd = os.getcwd()
    readme = os.path.join(wd, "README.txt")
    metadata, lines = genomepy.files.read_readme(readme)

    metadata["name"] = "my_cool_genome"
    lines = ["", "I wanted to do some regex, but regex is hard"]
    genomepy.files.write_readme(readme, metadata, lines)

    metadata2, lines2 = genomepy.files.read_readme(readme)
    assert metadata == metadata2
    assert lines == lines2

    os.unlink(readme)


def test_update_readme():
    wd = os.getcwd()
    readme = os.path.join(wd, "README.txt")

    # new data added
    updated_metadata = {"key": "value"}
    extra_lines = ["let's add some text to the lines section"]
    genomepy.files.update_readme(readme, updated_metadata, extra_lines)
    metadata, lines = genomepy.files.read_readme(readme)
    assert metadata["key"] == updated_metadata["key"]
    assert extra_lines[0] in lines

    # new data overwrites
    updated_metadata = {"key": "value2"}
    genomepy.files.update_readme(readme, updated_metadata)
    metadata, lines = genomepy.files.read_readme(readme)
    assert metadata["key"] == updated_metadata["key"]

    os.unlink(readme)


def test_generate_gap_bed(fname="tests/data/gap.fa", outname="tests/data/gap.bed"):
    tmp = NamedTemporaryFile().name
    genomepy.files.generate_gap_bed(fname, tmp)

    result = open(tmp).read()
    expected = open(outname).read()

    assert result == expected


def test_generate_fa_sizes(infa="tests/data/gap.fa"):
    tmp = NamedTemporaryFile().name
    genomepy.files.generate_fa_sizes(infa, tmp)

    result = open(tmp).read()

    assert result == "chr1\t28\nchr2\t45\nchr3\t15\n"


def test_tar_to_bigfile():
    fname = "tests/data/tar2.fa.tar.gz"
    outname = "tests/data/tar2.fa"
    genomepy.files.tar_to_bigfile(fname, outname)

    assert os.path.exists(outname)
    # tar2.fa is a copy of tar1.fa. Check if they are identical after untarring.
    assert filecmp.cmp(outname, "tests/data/tar1.fa")
    os.unlink(outname)


def test_gunzip_and_name(fname="tests/data/small_genome.fa.gz"):
    assert os.path.exists(fname)
    fname, gzip_file = genomepy.files.gunzip_and_name(fname)
    assert gzip_file and fname.endswith(".fa")
    assert os.path.exists(fname)
    assert not os.path.exists(fname + ".gz")


def test_gzip_and_name(fname="tests/data/small_genome.fa"):
    assert os.path.exists(fname)
    fname = genomepy.files.gzip_and_name(fname)
    assert fname.endswith(".gz")
    assert os.path.exists(fname)
    assert not os.path.exists(fname[:-3])

    fname, _ = genomepy.files.gunzip_and_name(fname)
    assert fname.endswith(".fa")
    assert os.path.exists(fname)
    assert not os.path.exists(fname + ".gz")


def test_bgzip_and_name(fname="tests/data/small_genome.fa"):
    assert os.path.exists(fname)
    fname = genomepy.files.bgzip_and_name(fname)
    assert fname.endswith(".gz")
    assert os.path.exists(fname)
    assert not os.path.exists(fname[:-3])

    with pytest.raises(sp.CalledProcessError):
        genomepy.files.bgzip_and_name("tests/data/nofile.fa")


def test_delete_extensions():
    fpath1 = "tests/data/weird_ext1.test123"
    fpath2 = "tests/data/weird_ext2.test123"
    for fpath in [fpath1, fpath2]:
        with open(fpath, "w") as f:
            f.write("asd\n")

    assert os.path.exists(fpath1)
    assert os.path.exists(fpath2)
    genomepy.files.delete_extensions("tests/data", ["test123"])
    assert not os.path.exists(fpath1)
    assert not os.path.exists(fpath2)


def test__open():
    # read/write regular file
    a1 = "tests/data/data.annotation.gtf"
    with genomepy.files._open(a1, "w") as gtf:
        gtf.write("regular file")

    with open(a1) as gtf:
        lines1 = gtf.readlines()
    assert lines1 == ["regular file"]

    with genomepy.annotation._open(a1) as gtf:
        lines2 = gtf.readlines()
    assert lines2 == lines1
    genomepy.files.rm_rf(a1)

    # read/write gzipped file
    a2 = "tests/data/data.annotation.gtf.gz"
    with genomepy.files._open(a2, "w") as gtf:
        gtf.write("gzipped file")

    with pytest.raises(UnicodeDecodeError):
        with open(a2) as gtf:
            gtf.read()

    with genomepy.annotation._open(a2) as gtf:
        lines = gtf.readlines()
    assert lines == ["gzipped file"]
    genomepy.files.rm_rf(a1)


def test_file_len():
    assert genomepy.files.file_len("tests/data/gap.fa") == 13


def test_get_file_info(fname="tests/data/small_genome.fa.gz"):
    ext, gz = genomepy.files.get_file_info(fname)
    assert ext == ".fa" and gz

    ext, gz = genomepy.files.get_file_info(fname[:-2] + "fai")
    assert ext == ".fai" and not gz


def test_glob_ext_files(file="tests/data/small_genome.fa"):
    assert file not in genomepy.files.glob_ext_files("tests/data")
    assert file + ".gz" in genomepy.files.glob_ext_files("tests/data")
    assert len(genomepy.files.glob_ext_files("tests/data", "fake_ext")) == 0
    assert len(genomepy.files.glob_ext_files("tests/data/regexp")) == 1


def test_is_genome_dir():
    # dir contains a fasta
    assert genomepy.files.is_genome_dir("tests/data/regexp")
    # dir does not contain a fasta
    assert not genomepy.files.is_genome_dir("tests/genome")


def test_filter_fasta(fname="tests/data/regexp/regexp.fa"):
    # no sequences left after filtering
    with pytest.raises(Exception), NamedTemporaryFile(suffix=".fa") as tmpfa:
        regex = "missing_chromosome"
        genomepy.files.filter_fasta(fname, regex=regex, outfa=tmpfa.name)

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
            keys = genomepy.files.filter_fasta(
                fname, regex=regex, invert_match=False, outfa=tmpfa.name
            ).keys()
            assert len(keys) == match, regex

        with NamedTemporaryFile(suffix=".fa") as tmpfa:
            keys = genomepy.files.filter_fasta(
                fname, regex=regex, invert_match=True, outfa=tmpfa.name
            ).keys()
            assert len(keys) == no_match, regex
