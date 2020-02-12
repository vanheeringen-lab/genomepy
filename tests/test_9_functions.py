import os
import pytest

from appdirs import user_config_dir
from shutil import rmtree
from tempfile import mkdtemp

import genomepy.functions


def test_manage_config():
    # will give a NameError if the config/config dir is missing
    rmtree(user_config_dir("genomepy"))
    for cmd in ["generate", "file", "show"]:
        genomepy.functions.manage_config(cmd)


def test_list_available_genomes(provider="ncbi"):
    g = genomepy.functions.list_available_genomes

    # any provider
    assert isinstance(next(g()), list)
    assert len(next(g())) == 3

    # specific provider
    assert next(g(provider))[0] == provider

    # wrong provider, throw exception
    with pytest.raises(Exception):
        next(g("not a provider"))


def test_list_available_providers():
    # NCBI, Ensembl, UCSC and direct URL (4 providers total)
    assert len(genomepy.functions.list_available_providers()) == 4


def test__is_genome_dir():
    # dir contains a fasta
    assert genomepy.functions._is_genome_dir("tests/data")
    # dir does not contain a fasta
    assert not genomepy.functions._is_genome_dir("tests/genome")


def test_list_installed_genomes():
    assert isinstance(genomepy.functions.list_installed_genomes(), list)


def test_search():
    # unrecognized provider/genome will cause an Exception or StopIteration respectively
    s = genomepy.functions.search
    assert isinstance(next(s("Xenopus Tropicalis", "NCBI")), list)
    assert isinstance(next(s("Xenopus Tropicalis")), list)


# install_genome is tested elsewhere


def test_generate_exports(genome="ASM14646v1", provider="NCBI"):
    tmp = mkdtemp()
    genomepy.install_genome(genome, provider, genome_dir=tmp)
    assert isinstance(genomepy.functions.generate_exports(), list)
    rmtree(tmp)


def test_generate_env():
    fname = "myenv.txt"
    genomepy.functions.generate_env(fname)

    config_dir = str(user_config_dir("genomepy"))
    path = os.path.join(config_dir, fname)
    assert os.path.exists(path)

    os.remove(path)


def test_glob_ext_files():
    gef = genomepy.functions.glob_ext_files
    assert "tests/data/small_genome.fa.gz" in gef("tests/data")


def test_genome():
    genomepy.functions.Genome("tests/data/gap.fa")


def test_manage_plugins():
    genomepy.functions.manage_plugins("list")
    genomepy.functions.manage_plugins("enable", ["star"])
    genomepy.functions.manage_plugins("disable", ["star"])
