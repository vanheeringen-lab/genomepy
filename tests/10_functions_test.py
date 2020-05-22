import genomepy
import pytest
import os
import shutil

from appdirs import user_config_dir
from platform import system

linux = system() == "Linux"
travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"


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
    assert genomepy.functions._is_genome_dir("tests/data")
    # dir does not contain a fasta
    assert not genomepy.functions._is_genome_dir("tests/genome")


def test_list_installed_genomes():
    assert isinstance(genomepy.functions.list_installed_genomes(os.getcwd()), list)

    gdir = os.path.join(os.getcwd(), "tests", "data")
    genomes = genomepy.functions.list_installed_genomes(gdir)
    assert genomes == ["regexp"]


@pytest.mark.skipif(not travis or not linux, reason="slow")
def test_install_genome():
    localname = "my_genome"
    genomepy.functions.install_genome(
        name="fr3",
        provider="UCSC",
        genomes_dir=None,
        localname=localname,
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
        genomes_dir, localname, localname + ".annotation.gtf.gz"
    )
    assert os.path.exists(annotation_file)

    readme = os.path.join(os.path.dirname(genome_file), "README.txt")
    with open(readme) as f:
        metadata = {}
        for line in f.readlines():
            vals = line.strip().split(":")
            metadata[vals[0].strip()] = (":".join(vals[1:])).strip()

    assert metadata["name"] == localname


@pytest.mark.skipif(
    not travis or not linux, reason="only works if a genome was installed"
)
def test_generate_exports():
    # already used, but we had to install a genome first to test it
    exports = genomepy.functions.generate_exports()
    assert isinstance(exports, list)
    # check if my_genome was installed in the last test
    assert any([x for x in exports if x.startswith("export MY_GENOME")])

    # add genome that throws a FastaIndexingError
    gd = genomepy.utils.get_genomes_dir(None, True)
    os.makedirs(os.path.join(gd, "testgenome"), exist_ok=True)
    path = os.path.join(gd, "testgenome", "testgenome.fa")
    with open(path, "w") as fa:
        fa.write("forbidden characters")
    exports = genomepy.functions.generate_exports()
    assert f"export TESTGENOME={path}" not in exports

    # add genome that works
    with open(path, "w") as fa:
        fa.write(">chr1\nallowed characters")
    exports = genomepy.functions.generate_exports()
    assert f"export TESTGENOME={path}" in exports

    shutil.rmtree(os.path.join(gd, "testgenome"))


@pytest.mark.skipif(
    not travis or not linux, reason="only works if a genome was installed"
)
def test_generate_env():
    # already used, but we had to install a genome first to test it
    config_dir = str(user_config_dir("genomepy"))
    path = os.path.join(config_dir, "exports.txt")
    if os.path.exists(path):
        os.unlink(path)
    assert not os.path.exists(path)

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
    search = genomepy.functions.search("xenopus tropicalis", "ensembl")
    metadata = next(search)
    assert isinstance(metadata, list)
    assert "xenopus_tropicalis_v9.1" in str(metadata[0])
    assert "ensembl" in str(metadata[1])
    assert "gca_000004195" in str(metadata[2])
    assert "8364" in str(metadata[4])


def test_accession_search():
    search = [row for row in genomepy.functions.search("GCA_000004195.3")]
    assert 3 == len(search)
    providers = [row[1] for row in search]
    assert b"Ensembl" in providers
    assert b"NCBI" in providers
    assert b"UCSC" in providers
