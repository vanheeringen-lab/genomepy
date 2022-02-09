import os
import subprocess as sp
from shutil import rmtree

import norns
import pandas as pd
import pytest
from appdirs import user_cache_dir
from diskcache import Cache

import genomepy

from . import linux, travis


def test_linting():
    out = sp.check_output(
        "chmod +x tests/format.sh; tests/format.sh lint",
        stderr=sp.STDOUT,  # send errors to out
        shell=True,
    )
    out = out.decode("utf-8").replace("\x1b[0m", "").replace("\nDone\n", "")
    if out != "":
        pytest.fail(
            f"Linting failed. Messages: \n\n{out}\n"
            "Run `tests/format.sh` to format the repo.",
            False,
        )


# caching.py


@pytest.mark.skipif(not travis, reason="it works locally all right")
def test_clean():
    genomepy.clean()
    my_cache_dir = os.path.join(user_cache_dir("genomepy"), genomepy.__version__)
    assert os.path.exists(my_cache_dir)  # dir exists
    assert not os.listdir(my_cache_dir)  # contains 0 pickles
    genomepy.clean()  # no errors when cache dir is empty


@pytest.mark.skipif(not travis, reason="it works locally all right")
def test_cache():
    # test caching of complex data types
    my_cache_dir = os.path.join(user_cache_dir("genomepy"), str(linux))
    os.makedirs(my_cache_dir, exist_ok=True)
    cache = Cache(directory=my_cache_dir, size_limit=1000000)
    test = ["a", "b", "c"]

    @cache.memoize(expire=10, tag="expensive_function")
    def expensive_function(data):
        return pd.DataFrame(data)

    # Get full function name https://github.com/grantjenks/python-diskcache/blob/master/diskcache/core.py
    def full_name(func):
        """Return full name of `func` by adding the module and function name."""
        return func.__module__ + "." + func.__qualname__

    # cache key tuple
    cache_key = (
        full_name(expensive_function),
        test,
        None,
    )

    # check that results before/after caching are identical
    expected = expensive_function(test)
    cached_data = cache.get(cache_key)
    assert cached_data.equals(expected), "Cached data does not match expected data"
    rmtree(my_cache_dir, ignore_errors=True)


# exceptions.py


def test_exceptions():
    with pytest.raises(genomepy.exceptions.GenomeDownloadError):
        raise genomepy.exceptions.GenomeDownloadError


# config.__init__.py


def test_config():
    config = norns.config("genomepy", default="config/default.yaml")
    assert len(config.keys()) == 3


def test_manage_config(capsys):
    # make a new config
    genomepy.config.manage_config("generate")
    captured = capsys.readouterr().out.strip()
    assert captured.startswith("Created config file")

    # check where it is found
    fname = os.path.expanduser("~/Library/Application Support/genomepy/genomepy.yaml")
    if linux:
        fname = os.path.expanduser("~/.config/genomepy/genomepy.yaml")
    genomepy.config.manage_config("file")
    captured = capsys.readouterr().out.strip()
    assert captured == fname

    # mess with the config
    with open(fname, "w") as f:
        print("bgzip: na", file=f)

    # show the mess
    genomepy.config.manage_config("show")
    captured = capsys.readouterr().out.strip()
    assert captured.startswith("bgzip: na")

    # make a new config
    genomepy.config.manage_config("generate")
    captured = capsys.readouterr().out.strip()
    assert captured.startswith("Created config file")

    # check if the mess was fixed
    genomepy.config.manage_config("show")
    captured = capsys.readouterr().out.strip()
    assert captured.startswith("bgzip: false")
