import genomepy
import norns
import os
import pytest
import subprocess as sp

from shutil import rmtree
from bucketcache import Bucket
from appdirs import user_cache_dir
from genomepy.__about__ import __version__

travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"


@pytest.mark.skipif(travis, reason="format before committing!")
def test_black_formatting():
    sp.check_call("black setup.py genomepy/ tests/", shell=True)


@pytest.mark.skipif(travis, reason="format before committing!")
def test_flake8_formatting():
    ret = sp.check_call("flake8 setup.py genomepy/ tests/", shell=True)
    assert ret == 0


def test_import():
    # __init__.py
    assert str(genomepy.search).startswith("<bound method")
    assert str(genomepy.Genome) == "<class 'genomepy.genome.Genome'>"
    assert str(genomepy.ProviderBase) == "<class 'genomepy.provider.ProviderBase'>"
    assert genomepy.__author__ == "Simon van Heeringen"


def test_exceptions():
    with pytest.raises(genomepy.exceptions.GenomeDownloadError):
        raise genomepy.exceptions.GenomeDownloadError


def test_config():
    config = norns.config("genomepy", default="cfg/default.yaml")
    assert len(config.keys()) == 3


@pytest.mark.skipif(not travis, reason="it works locally all right")
def test_cache(capsys):
    my_cache_dir = os.path.join(user_cache_dir("genomepy"), __version__)
    if os.path.exists(my_cache_dir):
        rmtree(my_cache_dir)
    os.makedirs(my_cache_dir)
    cached = Bucket(my_cache_dir, days=7)

    @cached
    def expensive_method():
        print("Method called.")

    @expensive_method.callback
    def expensive_method(callinfo):
        print("Cache used.")

    expensive_method()
    expensive_method()

    captured = capsys.readouterr().out.strip().split("\n")
    assert captured == ["Method called.", "Cache used."]
