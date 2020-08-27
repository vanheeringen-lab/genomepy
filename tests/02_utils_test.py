import filecmp
import genomepy
import pytest
import os
import shutil
import urllib.request
import subprocess as sp

from tempfile import NamedTemporaryFile
from platform import system

linux = system() == "Linux"


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
    tmpfa = NamedTemporaryFile(suffix=".fa").name

    # input is output error
    with pytest.raises(ValueError):
        genomepy.utils.filter_fasta(fname, fname)

    # output exists error
    with pytest.raises(FileExistsError):
        existing_file = "tests/data/gap.fa"
        genomepy.utils.filter_fasta(fname, existing_file)

    # no sequences left after filtering
    with pytest.raises(Exception):
        regex = "missing_chromosome"
        genomepy.utils.filter_fasta(fname, tmpfa, regex=regex)

    # function proper
    regexps = [
        ("Chr.*", 2, 15),
        ("Scaffold.*", 1, 16),
        ("scaffold_.*", 3, 14),
        (r"^\d+$", 4, 13),
        ("chr.*", 4, 13),
    ]
    for regex, match, no_match in regexps:
        fa = genomepy.utils.filter_fasta(fname, tmpfa, regex=regex, v=False, force=True)
        assert len(fa.keys()) == match
        fa = genomepy.utils.filter_fasta(fname, tmpfa, regex=regex, v=True, force=True)
        assert len(fa.keys()) == no_match


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


def test_cmd_ok():
    assert genomepy.utils.cmd_ok("STAR")
    assert genomepy.utils.cmd_ok("bwa")
    assert not genomepy.utils.cmd_ok("missing_cmd")


def test_run_index_cmd(capsys, name="tests", good_cmd="ls", bad_cmd="bad_cmd"):
    genomepy.utils.run_index_cmd(name=name, cmd=good_cmd)

    # bad_command not found error
    genomepy.utils.run_index_cmd(name=name, cmd=bad_cmd)
    captured = capsys.readouterr().err.strip()
    if linux:
        assert captured.endswith(f"{bad_cmd}: not found")
    else:
        # thanks for being different mac. That took me 30 minutes of my life...
        assert captured.endswith(f"{bad_cmd}: command not found")


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

    # URL input
    url = "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_genome.fa.gz"
    result = genomepy.utils.get_localname(name=url)
    assert result == name


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


def test_is_number():
    assert genomepy.utils.is_number("1234")
    assert genomepy.utils.is_number(1234)
    assert not genomepy.utils.is_number("abcd")


def test_retry(capsys):
    # runs passed function
    txt = "hello world"
    genomepy.utils.retry(print, 1, txt)
    captured = capsys.readouterr().out.strip()
    assert captured == txt

    # handles URLErrors
    def _offline_func():
        raise urllib.request.URLError("this function is offline")

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


def test_match_contigs():
    sizes_file = "tests/data/small_genome.fa.sizes"
    gtf_file = "tests/data/small_genome.annotation.gtf"

    # generate sizes file
    with open(sizes_file, "w") as f:
        f.write("""chr2\t1234""")

    # generate empty gtf file
    with open(gtf_file, "w") as f:
        f.write("# skip this line\n")
    match = genomepy.utils.match_contigs(gtf_file, sizes_file)
    assert match is None

    # generate matching gtf file
    with open(gtf_file, "w") as f:
        f.write("# skip this line\n")
        f.write("""chr2\tvanHeeringen-lab""")
    match = genomepy.utils.match_contigs(gtf_file, sizes_file)
    assert match

    # generate mismatching gtf file
    with open(gtf_file, "w") as f:
        f.write("# skip this line\n")
        f.write("""2\tvanHeeringen-lab""")
    match = genomepy.utils.match_contigs(gtf_file, sizes_file)
    assert match is False


def test_contig_pos():
    gtf_file = "tests/data/small_genome.annotation.gtf"
    genome_file = "tests/data/contig_pos.fa"

    # generate gtf file
    with open(gtf_file, "w") as f:
        f.write("chr2\tvanHeeringen-lab\n" "chr3\tvanHeeringen-lab\n")

    # generate non-matching genome file
    with open(genome_file, "w") as f:
        f.write(">no matches\n")
    element_pos = genomepy.utils.contig_pos(gtf_file, genome_file)
    assert element_pos == -1

    # generate matching genome file
    with open(genome_file, "w") as f:
        f.write(">chr3 this matches\n")
    element_pos = genomepy.utils.contig_pos(gtf_file, genome_file)
    assert element_pos == 0

    # generate different matching genome file
    with open(genome_file, "w") as f:
        f.write(">this matches chr2\n")
    element_pos = genomepy.utils.contig_pos(gtf_file, genome_file)
    assert element_pos == 2

    os.unlink(genome_file)


def test_contig_conversion():
    genome_file = "tests/data/contig_pos.fa"

    # generate genome file
    with open(genome_file, "w") as f:
        f.write(">this matches chr2\n" ">that matches chr3\n")
    ids, duplicate_contigs = genomepy.utils.contig_conversion(genome_file, 2)
    assert not duplicate_contigs
    assert ids["chr2"] == "this"
    assert ids["chr3"] == "that"

    # assume "matches" is the contig name
    ids, duplicate_contigs = genomepy.utils.contig_conversion(genome_file, 1)
    assert "matches" in duplicate_contigs
    assert ids["matches"] == "that"

    os.unlink(genome_file)


def test_sanitize_gtf():
    old_gtf = "tests/data/old_gtf.gtf"
    new_gtf = "tests/data/new_gtf.gtf"
    ids = {"1": "chr1", "2": "chr2"}  # replace 1,2 with chr1, chr2

    with open(old_gtf, "w") as gtf:
        gtf.write("1\tgene1\n" "2\tgene2\n" "3\tgene3\n")
    result = open(old_gtf).read()
    assert result == "1\tgene1\n2\tgene2\n3\tgene3\n"

    missing_contigs = genomepy.utils.sanitize_gtf(old_gtf, new_gtf, ids)
    result = open(new_gtf).read()
    assert result == "chr1\tgene1\nchr2\tgene2\n3\tgene3\n"
    assert missing_contigs == ["3"]

    os.unlink(old_gtf)
    os.unlink(new_gtf)


def test_sanitize_annotation(genome="tests/data/small_genome.fa.gz"):
    class TestGenome:
        def __init__(self, filename=genome):
            self.filename = filename
            self.genome_dir = os.path.dirname(genome)
            self.annotation_gtf_file = genome[:-5] + "annotation.gtf.gz"
            self.annotation_bed_file = genome[:-5] + "annotation.bed.gz"
            self.sizes_file = genome[:-2] + "sizes"
            self.readme_file = os.path.join(self.genome_dir, "README.txt")

    # generate sizes file (already tested)
    sp.check_call(f"gunzip -f {genome}", shell=True)
    sizes_file = genome[:-2] + "sizes"
    genomepy.utils.generate_fa_sizes(genome[:-3], sizes_file)
    sp.check_call(f"bgzip -f {genome[:-3]}", shell=True)

    # generate gtf file
    gtf_file = genome[:-5] + "annotation.gtf"
    with open(gtf_file, "w") as f:
        f.write("# skip this line\n")
        f.write(
            """chr2\tvanHeeringen-lab\tgene\t2\t22\t.\t+\t.\tgene_id "vH-1"; transcript_id "vH-1.1";\n"""
        )
    sp.check_call(f"gzip -f {gtf_file}", shell=True)

    # generate bed file
    bed_file = gtf_file.replace("gtf", "bed.gz")
    sp.check_call(f"touch {bed_file}", shell=True)

    # function proper
    g = TestGenome(genome)
    genomepy.utils.sanitize_annotation(g)

    sp.check_call(f"gunzip -f {gtf_file}.gz", shell=True)
    result = open(gtf_file).read()
    assert result.startswith("# skip this line\nchr2\tvanHeeringen-lab")

    # cleanup
    os.unlink(sizes_file)
    os.unlink(gtf_file)
    os.unlink(bed_file)
