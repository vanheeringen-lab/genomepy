import pytest
import os.path
from appdirs import user_config_dir

import genomepy.functions


def test_manage_config():
    # will give a NameError if the config/config dir is missing
    genomepy.functions.manage_config("file")
    genomepy.functions.manage_config("show")


def test_list_available_genomes(provider="ncbi"):
    g = genomepy.functions.list_available_genomes(provider)
    for row in g:
        assert ("\t".join(row)).startswith(provider)
        break

    with pytest.raises(Exception):
        g = genomepy.functions.list_available_genomes('not a provider')
        for row in g:
            print("\t".join(row))
            break


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
    assert isinstance(next(genomepy.functions.search("Xenopus Tropicalis", "NCBI")), list)


# skipping several large/vague functions


def test_generate_exports():
    assert isinstance(genomepy.functions.generate_exports(), list)


# def test_generate_env():
#     genomepy.functions.generate_env("tests/data/gap.fa")
#
#     config_dir = str(user_config_dir("genomepy"))
#     path = os.path.join(config_dir, "exports.txt")
#     assert os.path.exists(path)
#
#     os.remove(path)


def test_glob_ext_files():
    assert 'tests/data/small_genome.fa.gz' in genomepy.functions.glob_ext_files("tests/data")


def test_genome():
    genomepy.functions.Genome("tests/data/gap.fa")


def test_manage_plugins():
    genomepy.functions.manage_plugins("list")
