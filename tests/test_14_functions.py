import os

import pytest
from appdirs import user_config_dir

import genomepy
import genomepy.utils
from tests import linux, travis


def test_clean():
    # test moved to 01_tests to prevent errors in parallel tests
    pass


def test_manage_config(capsys):
    # make a new config
    genomepy.functions.manage_config("generate")
    captured = capsys.readouterr().out.strip()
    assert captured.startswith("Created config file")

    # check where it is found
    fname = os.path.expanduser("~/Library/Application Support/genomepy/genomepy.yaml")
    if linux:
        fname = os.path.expanduser("~/.config/genomepy/genomepy.yaml")
    genomepy.functions.manage_config("file")
    captured = capsys.readouterr().out.strip()
    assert captured == fname

    # mess with the config
    with open(fname, "w") as f:
        print("bgzip: na", file=f)

    # show the mess
    genomepy.functions.manage_config("show")
    captured = capsys.readouterr().out.strip()
    assert captured.startswith("bgzip: na")

    # make a new config
    genomepy.functions.manage_config("generate")
    captured = capsys.readouterr().out.strip()
    assert captured.startswith("Created config file")

    # check if the mess was fixed
    genomepy.functions.manage_config("show")
    captured = capsys.readouterr().out.strip()
    assert captured.startswith("bgzip: false")


def test_list_available_genomes():
    g = genomepy.functions.list_available_genomes("Ensembl")
    metadata = next(g)
    assert isinstance(metadata, list)
    assert metadata[0:2] == ["athCun1", "Ensembl"]


def test_list_installed_genomes():
    assert isinstance(genomepy.functions.list_installed_genomes(os.getcwd()), list)

    gdir = os.path.join(os.getcwd(), "tests", "data")
    genomes = genomepy.functions.list_installed_genomes(gdir)
    assert set(genomes) == {"regexp", "sacCer3"}  # OSX likes to sort differently

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

    # find the first provider (Ensembl)
    provider = None
    p = genomepy.functions._lazy_provider_selection(name, provider)
    assert "ensembl" in str(p)

    # cant find genome anywhere
    name = "not_a_genome"
    with pytest.raises(genomepy.GenomeDownloadError):
        genomepy.functions._lazy_provider_selection(name, provider)


def test__provider_selection():
    # specified provider
    name = "Xenopus_tropicalis_v9.1"
    localname = "test_genome"
    genomes_dir = os.getcwd()
    provider = "NCBI"
    p = genomepy.functions._provider_selection(name, localname, genomes_dir, provider)
    assert "ncbi" in str(p)

    # provider from readme
    readme = os.path.join(genomes_dir, localname, "README.txt")
    os.makedirs(os.path.dirname(readme), exist_ok=True)
    with open(readme, "w") as r:
        r.write("provider: NCBI")
    provider = None
    p = genomepy.functions._provider_selection(name, localname, genomes_dir, provider)
    assert "ncbi" in str(p)
    genomepy.utils.rm_rf(os.path.dirname(readme))

    # lazy provider
    p = genomepy.functions._provider_selection(name, localname, genomes_dir, provider)
    assert "ensembl" in str(p)


def test__filter_genome():
    pass  # TODO


@pytest.mark.skipif(not travis, reason="slow")
def test_install_genome():
    localname = "my_genome"
    genomepy.functions.install_genome(
        name="dm3",
        provider="UCSC",
        genomes_dir=None,
        localname=localname,
        regex="R",
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
    assert "chr2R" in sizes
    assert "chr2L" not in sizes


# already used, but we had to install a genome first to test it
@pytest.mark.skipif(not travis, reason="a genome must be installed")
def test_generate_exports():
    exports = genomepy.functions.generate_exports()
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
    exports = genomepy.functions.generate_exports()
    assert f"export TESTGENOME={path}" not in exports

    # add genome that works
    with open(path, "w") as fa:
        fa.write(">chr1\nallowed characters")
    genomepy.Genome(name, gd)  # create index
    exports = genomepy.functions.generate_exports()
    assert f"export TESTGENOME={path}" in exports

    genomepy.utils.rm_rf(os.path.join(gd, "testgenome"))


# already used, but we had to install a genome first to test it
@pytest.mark.skipif(not travis, reason="a genome must be installed")
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


def test_manage_plugins(capsys):
    genomepy.functions.manage_plugins("enable", ["blacklist"])
    genomepy.functions.manage_plugins("list")
    captured = capsys.readouterr().out.strip().split("\n")
    assert captured[2].startswith("blacklist")
    assert captured[2].endswith("*")

    genomepy.functions.manage_plugins("disable", ["blacklist"])
    genomepy.functions.manage_plugins("list")
    captured = capsys.readouterr().out.strip().split("\n")
    assert captured[2].startswith("blacklist")
    assert not captured[2].endswith("*")

    with pytest.raises(ValueError):
        genomepy.functions.manage_plugins("blurp")
