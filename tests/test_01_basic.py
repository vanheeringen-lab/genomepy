from appdirs import user_cache_dir
from bucketcache import Bucket
import os
import subprocess as sp
from shutil import rmtree

import norns
import pytest

import genomepy
from . import travis


def test_flake_formatting():
    """remove unused stuff (ignores lines marked with '# noqa')"""
    try:
        sp.check_output(
            "autoflake -r "
            f'{"--check " if travis else "--in-place "}'
            "--remove-all-unused-imports "
            "--remove-duplicate-keys "
            "--remove-unused-variables "
            "setup.py genomepy/ tests/",
            stderr=sp.STDOUT,
            shell=True,
        )
    except sp.CalledProcessError as e:
        msg = e.output.decode("utf-8")
        msg = msg.replace("No issues detected!", "")
        pytest.fail(msg, False)


def test_isort_formatting():
    """sort imports"""
    try:
        sp.check_output(
            "isort "
            f'{"--check " if travis else "--overwrite-in-place "}'
            "--profile black "
            "--conda-env environment.yml "
            "setup.py genomepy/ tests/",
            stderr=sp.STDOUT,
            shell=True,
        )
    except sp.CalledProcessError as e:
        msg = e.output.decode("utf-8")
        pytest.fail(msg, False)


def test_black_formatting():
    try:
        sp.check_output(
            f"black {'--check ' if travis else ''} setup.py genomepy/ tests/",
            stderr=sp.STDOUT,
            shell=True,
        )
    except sp.CalledProcessError as e:
        msg = e.output.decode("utf-8")
        msg = msg.split("\n")[:-3]
        msg = "\n".join(["Black output:"] + msg)
        pytest.fail(msg, False)


def test_flake8_linting():
    try:
        sp.check_output(
            "flake8 setup.py genomepy/ tests/", stderr=sp.STDOUT, shell=True
        )
    except sp.CalledProcessError as e:
        msg = e.output.decode("utf-8")
        msg = "Flake8 output:\n" + msg
        pytest.fail(msg, False)


@pytest.mark.skipif(not travis, reason="it works locally all right")
def test_clean():
    # test moved from 10 to prevent errors in parallel tests
    genomepy.clean()

    my_cache_dir = os.path.join(user_cache_dir("genomepy"), genomepy.__version__)
    assert os.path.exists(my_cache_dir)  # dir exists
    assert not os.listdir(my_cache_dir)  # contains 0 pickles
    genomepy.clean()  # no errors when cache dir is empty


def test_import():
    # __init__.py
    assert str(genomepy.Genome) == "<class 'genomepy.genome.Genome'>"
    assert str(genomepy.Provider) == "<class 'genomepy.provider.Provider'>"
    assert "Simon van Heeringen" in genomepy.__author__


def test_exceptions():
    with pytest.raises(genomepy.exceptions.GenomeDownloadError):
        raise genomepy.exceptions.GenomeDownloadError


def test_config():
    config = norns.config("genomepy", default="cfg/default.yaml")
    assert len(config.keys()) == 3


@pytest.mark.skipif(not travis, reason="it works locally all right")
def test_cache(capsys):
    my_cache_dir = os.path.join(user_cache_dir("genomepy"), genomepy.__version__)
    if os.path.exists(my_cache_dir):
        rmtree(my_cache_dir)
    os.makedirs(my_cache_dir)
    cache = Bucket(my_cache_dir, days=7)

    @cache
    def expensive_method():
        print("Method called.")

    @expensive_method.callback
    def expensive_method(_):
        print("Cache used.")

    expensive_method()
    expensive_method()

    captured = capsys.readouterr().out.strip().split("\n")
    assert captured == ["Method called.", "Cache used."]
