import os

import pytest
from appdirs import user_config_dir

import genomepy
import genomepy.utils
from tests import linux, travis


def test_head_annotations(caplog, capsys):
    genomepy.functions.head_annotations("ASM14646v1", provider="ncbi", n=1)
    captured = capsys.readouterr().out.strip()

    assert "NCBI" in caplog.text
    assert 'gene_name "Eint_010010";' in captured


def test_list_available_genomes():
    g = genomepy.functions.list_available_genomes("Ensembl")
    metadata = next(g)
    assert isinstance(metadata, list)
    assert metadata[0:2] == ["ASM394721v1", "Ensembl"]


def test_list_installed_genomes():
    assert isinstance(genomepy.functions.list_installed_genomes(os.getcwd()), list)

    gdir = os.path.join(os.getcwd(), "tests", "data")
    genomes = genomepy.functions.list_installed_genomes(gdir)
    assert set(genomes) == {
        "regexp",
        "sacCer3",
        "sanitize",
    }  # OSX likes to sort differently

    empty_list = genomepy.functions.list_installed_genomes("./thisdirdoesnotexist")
    assert empty_list == []


def test__lazy_provider_selection():
    # Xenopus_tropicalis_v9.1 can be found on both Ensembl and NCBI.
    # Ensembl is first in lazy selection.

    # find genome in specified provider (NCBI)
    name = "Xenopus_tropicalis_v9.1"
    provider = "NCBI"
    p = genomepy.functions._lazy_provider_selection(name, provider)
    assert "ncbi" in str(p)

    # GENCODE's FTP does not work on Travis-Linux
    if not (travis and linux):
        # find the first provider (Ensembl)
        provider = None
        p = genomepy.functions._lazy_provider_selection(name, provider)
        assert "ensembl" in str(p)

        # cant find genome anywhere
        name = "not_a_genome"
        with pytest.raises(genomepy.exceptions.GenomeDownloadError):
            genomepy.functions._lazy_provider_selection(name, provider)


def test__provider_selection():
    # specified provider
    name = "Xenopus_tropicalis_v9.1"
    localname = "test_genome"
    genomes_dir = os.getcwd()
    provider = "NCBI"
    p = genomepy.functions._provider_selection(name, localname, genomes_dir, provider)
    assert "ncbi" in str(p)

    # lazy provider
    p = genomepy.functions._provider_selection(name, localname, genomes_dir)
    assert "ensembl" in str(p)


def test__provider_selection_readme():
    name = "Xenopus_tropicalis_v9.1"
    localname = "test_genome"
    genomes_dir = os.getcwd()
    provider = "NCBI"
    # provider from readme
    readme = os.path.join(genomes_dir, localname, "README.txt")
    os.makedirs(os.path.dirname(readme), exist_ok=True)
    with open(readme, "w") as r:
        r.write("provider: NCBI")
    provider = None
    p = genomepy.functions._provider_selection(name, localname, genomes_dir, provider)
    assert "ncbi" in str(p)
    genomepy.utils.rm_rf(os.path.dirname(readme))


def test__get_fasta_regex_func():
    # filter alt regions (default)
    func = genomepy.functions._get_fasta_regex_func(regex=None, keep_alt=False)
    assert func("alt1") is False
    assert func("chr1") is True
    assert func("ALT1") is False  # case insensitive

    # filter user specified regex
    func = genomepy.functions._get_fasta_regex_func(
        regex="chr", invert_match=False, keep_alt=True
    )
    assert func("chr1") is True
    assert func("alt1") is False
    assert func("something_else") is False

    # filter user specified regex (inverted)
    func = genomepy.functions._get_fasta_regex_func(
        regex="chr", invert_match=True, keep_alt=True
    )
    assert func("chr1") is False
    assert func("alt1") is True
    assert func("something_else") is True

    # filter both
    func = genomepy.functions._get_fasta_regex_func(
        regex="chr", invert_match=True, keep_alt=False
    )
    assert func("chr1") is False
    assert func("alt1") is False
    assert func("something_else") is True


def test_install_genome():
    localname = "my_genome"
    genomepy.functions.install_genome(
        name="tests/data/sacCer3/sacCer3.fa",
        provider="Local",
        genomes_dir=None,
        localname=localname,
        regex="chrIV",
        annotation=True,
        force=True,
    )

    genomes_dir = genomepy.functions.get_genomes_dir(None, False)
    genome_file = os.path.join(genomes_dir, localname, localname + ".fa")
    assert os.path.exists(genome_file)
    sizes_file = os.path.join(genomes_dir, localname, localname + ".fa.sizes")
    assert os.path.exists(sizes_file)
    gaps_file = os.path.join(genomes_dir, localname, localname + ".gaps.bed")
    assert os.path.exists(gaps_file)
    annotation_file = os.path.join(
        genomes_dir, localname, localname + ".annotation.gtf"
    )
    assert os.path.exists(annotation_file)

    # regex test:
    sizes = genomepy.Genome(localname).sizes.keys()
    assert "chrIV" in sizes


# already used, but we had to install a genome first to test it
def test_generate_exports():
    exports = genomepy.functions._generate_exports()
    assert isinstance(exports, list)
    # check if my_genome was installed in the last test
    assert any([x for x in exports if x.startswith("export MY_GENOME")])

    # add genome that throws an IndexNotFoundError
    gd = genomepy.utils.get_genomes_dir(None, True)
    name = "testgenome"
    os.makedirs(os.path.join(gd, name), exist_ok=True)
    path = os.path.join(gd, name, f"{name}.fa")
    with open(path, "w") as fa:
        fa.write("genome without index")
    exports = genomepy.functions._generate_exports()
    assert f"export TESTGENOME={path}" not in exports

    # add genome that works
    with open(path, "w") as fa:
        fa.write(">chr1\nallowed characters")
    genomepy.Genome(name, gd)  # create index
    exports = genomepy.functions._generate_exports()
    assert f"export TESTGENOME={path}" in exports

    genomepy.utils.rm_rf(os.path.join(gd, "testgenome"))


# already used, but we had to install a genome first to test it
def test_generate_env():
    config_dir = str(user_config_dir("genomepy"))
    path = os.path.join(config_dir, "exports.txt")

    # give file path
    my_path = "~/exports.txt"
    genomepy.functions.generate_env(my_path)
    assert os.path.exists(os.path.expanduser(my_path))
    os.unlink(os.path.expanduser(my_path))

    # give file name
    my_file = os.path.join(config_dir, "my_exports.txt")
    genomepy.functions.generate_env("my_exports.txt")
    assert os.path.exists(my_file)
    os.unlink(os.path.expanduser(my_file))

    # give nothing
    if os.path.exists(path):
        os.unlink(path)
    genomepy.functions.generate_env()
    assert os.path.exists(path)

    with open(path) as f:
        exports = []
        for line in f.readlines():
            vals = line.strip()
            exports.append(vals)
    assert any([x for x in exports if x.startswith("export MY_GENOME")])
    os.unlink(path)


def test__delete_extensions():
    fpath1 = "tests/data/empty/weird_ext1.test123"
    fpath2 = "tests/data/empty/weird_ext2.test123"
    for fpath in [fpath1, fpath2]:
        with open(fpath, "w") as f:
            f.write("asd\n")

    assert os.path.exists(fpath1)
    assert os.path.exists(fpath2)
    genomepy.functions._delete_extensions("tests/data/empty", ["test123"])
    assert not os.path.exists(fpath1)
    assert not os.path.exists(fpath2)


def test__is_genome_dir():
    # dir contains a fasta
    assert genomepy.functions._is_genome_dir("tests/data/regexp")
    # dir does not contain a fasta
    assert not genomepy.functions._is_genome_dir("tests/genome")
