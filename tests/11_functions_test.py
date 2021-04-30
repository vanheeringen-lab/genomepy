import genomepy
import genomepy.utils
import genomepy.argparse_support
import pytest
import os
import argparse
from tempfile import TemporaryDirectory

from appdirs import user_config_dir
from platform import system

linux = system() == "Linux"
travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"


def test_clean():
    # test moved to 01_tests to prevent errors in parallel tests
    pass

    # my_cache_dir = os.path.join(
    #     user_cache_dir("genomepy"), genomepy.__about__.__version__
    # )
    #
    # genomepy.provider.ProviderBase.create("UCSC")  # pickles UCSC genomes
    # assert os.path.exists(my_cache_dir)  # dir exists
    # assert os.listdir(my_cache_dir)  # contains >=1 pickle(s)
    #
    # genomepy.clean()
    # assert os.path.exists(my_cache_dir)  # dir exists
    # assert not os.listdir(my_cache_dir)  # contains 0 pickles
    #
    # genomepy.clean()  # no errors when cache dir is empty


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


def test_online_providers():
    ops = genomepy.functions.online_providers()
    providers = [p.name for p in ops]
    assert len(providers) == 4
    assert providers[-1] == "URL"

    ops = genomepy.functions.online_providers("Ensembl")
    providers = [p.name for p in ops]
    assert len(providers) == 1
    assert providers[0] == "Ensembl"


def test_list_available_genomes():
    g = genomepy.functions.list_available_genomes("Ensembl")
    metadata = next(g)
    assert isinstance(metadata, list)
    assert metadata[0] == "Ensembl"

    # use a loop in case more genomes are added
    for genome in g:
        if genome[1] == "JCVI-ESG2-1.0":
            assert genome == [
                "Ensembl",
                "JCVI-ESG2-1.0",
                "GCA_000208925.2",
                "Entamoeba histolytica",
                "294381",
                "AmoebaDB_1.6",
            ]
            break

    g = genomepy.functions.list_available_genomes()
    for genome in g:
        if genome[1] == "ENA_1":
            assert genome == [
                "Ensembl",
                "ENA_1",
                "na",
                "Albugo laibachii",
                "890382",
                "2011-08-ENA",
            ]
            break


def test__is_genome_dir():
    # dir contains a fasta
    assert genomepy.functions._is_genome_dir("tests/data/regexp")
    # dir does not contain a fasta
    assert not genomepy.functions._is_genome_dir("tests/genome")


def test_list_installed_genomes():
    assert isinstance(genomepy.functions.list_installed_genomes(os.getcwd()), list)

    gdir = os.path.join(os.getcwd(), "tests", "data")
    genomes = genomepy.functions.list_installed_genomes(gdir)
    assert genomes == ["regexp"]

    empty_list = genomepy.functions.list_installed_genomes("./thisdirdoesnotexist")
    assert empty_list == []


def test__lazy_provider_selection():
    # Xenopus_tropicalis_v9.1 can be found on both Ensembl and NCBI.
    # Ensembl is first in lazy selection.

    # find genome in specified provider (NCBI)
    name = "Xenopus_tropicalis_v9.1"
    provider = "NCBI"
    p = genomepy.functions._lazy_provider_selection(name, provider)
    assert "NcbiProvider" in str(p)

    # find the first provider (Ensembl)
    provider = None
    p = genomepy.functions._lazy_provider_selection(name, provider)
    assert "EnsemblProvider" in str(p)

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
    assert "NcbiProvider" in str(p)

    # provider from readme
    readme = os.path.join(genomes_dir, localname, "README.txt")
    os.makedirs(os.path.dirname(readme), exist_ok=True)
    with open(readme, "w") as r:
        r.write("provider: NCBI")
    provider = None
    p = genomepy.functions._provider_selection(name, localname, genomes_dir, provider)
    assert "NcbiProvider" in str(p)
    genomepy.utils.rm_rf(os.path.dirname(readme))

    # lazy provider
    p = genomepy.functions._provider_selection(name, localname, genomes_dir, provider)
    assert "EnsemblProvider" in str(p)


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


@pytest.mark.skipif(not travis, reason="a genome must be installed")
def test_generate_exports():
    # already used, but we had to install a genome first to test it
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


@pytest.mark.skipif(not travis, reason="a genome must be installed")
def test_generate_env():
    # already used, but we had to install a genome first to test it
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


def test_list_available_providers():
    # NCBI, Ensembl, UCSC and direct URL (4 providers total)
    assert len(genomepy.functions.list_available_providers()) == 4


def test_search():
    # unrecognized provider/genome will cause an exception or stopiteration respectively
    # case insensitive description search
    search = genomepy.functions.search("xEnOpUs TrOpIcAlIs", "ensembl")
    metadata = next(search)

    # case insensitive assembly name search
    search = genomepy.functions.search("XeNoPuS_tRoPiCaLiS_v9.1", "ensembl")
    metadata2 = next(search)

    assert metadata == metadata2
    assert isinstance(metadata, list)
    assert "Xenopus_tropicalis_v9.1" in str(metadata[0])
    assert "Ensembl" in str(metadata[1])
    assert "GCA_000004195" in str(metadata[2])
    assert "8364" in str(metadata[4])


def test_accession_search():
    search = [row for row in genomepy.functions.search("GCA_000004195.3")]
    assert 3 == len(search)
    providers = [row[1] for row in search]
    assert b"Ensembl" in providers
    assert b"NCBI" in providers
    assert b"UCSC" in providers


def test_argparse_plugin():
    action = genomepy.argparse_support.parse_genome

    with TemporaryDirectory() as tmpdir:
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "-g", dest="genome", action=action(auto_install=True, genomes_dir=tmpdir)
        )
        args = parser.parse_args(
            ["-g", "ASM2732v1"],
        )
        assert isinstance(args.genome, genomepy.Genome)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        parser = argparse.ArgumentParser()
        parser.add_argument("-g", dest="genome", action=action())
        _ = parser.parse_args(["-g", "non_existing"])
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1
