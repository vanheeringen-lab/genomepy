import genomepy
import norns
import os
import pytest
import subprocess as sp

travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"


@pytest.mark.skipif(travis, reason="format before committing!")
def test_black_formatting():
    sp.check_call("black setup.py genomepy/ tests/", shell=True)


def test_import():
    # __init__.py
    assert str(genomepy.search).startswith("<function search at")
    assert str(genomepy.Genome) == "<class 'genomepy.genome.Genome'>"
    assert str(genomepy.ProviderBase) == "<class 'genomepy.provider.ProviderBase'>"
    assert genomepy.__author__ == "Simon van Heeringen"


def test_exceptions():
    with pytest.raises(genomepy.exceptions.GenomeDownloadError):
        raise genomepy.exceptions.GenomeDownloadError


def test_config():
    config = norns.config("genomepy", default="cfg/default.yaml")
    assert len(config.keys()) == 3
