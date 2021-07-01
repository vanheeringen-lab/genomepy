import os
import subprocess as sp
from shutil import rmtree

import norns
import pytest
from appdirs import user_cache_dir
from bucketcache import Bucket

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
def test_cache(capsys):
    my_cache_dir = os.path.join(user_cache_dir("genomepy"), str(linux))
    os.makedirs(my_cache_dir, exist_ok=True)
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
    rmtree(my_cache_dir, ignore_errors=True)
    assert captured == ["Method called.", "Cache used."]


# exceptions.py


def test_exceptions():
    with pytest.raises(genomepy.exceptions.GenomeDownloadError):
        raise genomepy.exceptions.GenomeDownloadError


# config.__init__.py


def test_config():
    config = norns.config("genomepy", default="config/default.yaml")
    assert len(config.keys()) == 3


def test_manage_config(caplog):
    # make a new config
    genomepy.config.manage_config("generate")
    assert "Created config file" in caplog.text

    # check where it is found
    fname = os.path.expanduser("~/Library/Application Support/genomepy/genomepy.yaml")
    if linux:
        fname = os.path.expanduser("~/.config/genomepy/genomepy.yaml")
    genomepy.config.manage_config("file")
    assert fname in caplog.text

    # mess with the config
    with open(fname, "w") as f:
        print("bgzip: na", file=f)

    # # show the mess
    # genomepy.config.manage_config("show")
    # assert "bgzip: na" in caplog.text
    # assert "bgzip: false" not in caplog.text
    #
    # # check if the mess was fixed
    # genomepy.config.manage_config("show")
    # assert "bgzip: false" in caplog.text
    # TODO: fix
