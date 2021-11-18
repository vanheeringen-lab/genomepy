"""Utility functions"""
import itertools
import os
import re
import shutil
import subprocess as sp
import sys
import time
from typing import Any, Optional
from urllib.request import urlopen

from loguru import logger

from genomepy.config import config


def cleanpath(path):
    """Expand any path input to a literal path output"""
    return os.path.abspath(  # expand relative paths ('./' and '../')
        os.path.expanduser(  # expand '~'
            os.path.expandvars(path)  # expand '$VARIABLES'
        )
    )


def mkdir_p(path):
    """'mkdir -p' in Python"""
    path = cleanpath(path)
    os.makedirs(path, exist_ok=True)


def rm_rf(path):
    """'rm -rf' in Python"""
    path = cleanpath(path)
    if os.path.isfile(path):
        try:
            os.unlink(path)
        except OSError:  # in case of NTFS related issues
            pass
    elif os.path.isdir(path):
        shutil.rmtree(path, ignore_errors=True)


def get_genomes_dir(genomes_dir: str = None, check_exist: Optional[bool] = True) -> str:
    """import genomes_dir if none is given, and check validity"""
    if not genomes_dir:
        # backwards compatibility for "genome_dir" (this fixes issue #87)
        genomes_dir = config.get("genomes_dir", config.get("genome_dir", None))
    if not genomes_dir:
        raise FileNotFoundError("Please provide or configure a genomes_dir")

    genomes_dir = cleanpath(genomes_dir)
    if not os.path.exists(genomes_dir) and check_exist:
        raise FileNotFoundError(f"Genomes_dir {genomes_dir} does not exist!")

    return genomes_dir


def cmd_ok(cmd) -> bool:
    """Returns True if cmd can be run."""
    try:
        sp.check_call(cmd, stderr=sp.PIPE, stdout=sp.PIPE)
    except sp.CalledProcessError:
        # bwa gives return code of 1 with no argument
        pass
    except FileNotFoundError:
        logger.error(f"{cmd} not found, skipping")
        return False
    return True


def run_index_cmd(name, cmd):
    """Run command, show errors if the returncode is non-zero."""
    logger.info(f"Creating {name} index...")
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)

    # show a spinner while the command is running
    spinner = itertools.cycle(["-", "\\", "|", "/"])
    while p.poll() is None:
        sys.stdout.write("\r" + next(spinner))
        time.sleep(0.15)
        sys.stdout.flush()
    sys.stdout.write("\b")  # clear the spinner

    stdout, stderr = p.communicate()
    if p.returncode != 0:
        logger.error(
            f"Indexing failed\n"
            f"stdout: {stdout.decode('utf8')}\n"
            f"stderr: {stderr.decode('utf8')}"
        )


def get_genomename(name):
    """return the name of the genome without path or extensions"""
    name = os.path.basename(name)  # remove path
    name = re.split(r"\.fa|\.fna|\.annot", name)[0]  # remove extensions
    return name


def get_remotename(name):
    """try to get the name from the url"""
    name = name.split("/")[-1]  # remove path
    name = name.replace(".gz", "")  # remove .gz
    name = os.path.splitext(name)[0]  # remove .fa/.fna/.fasta etc
    # remove unwanted substrings from the name (ex: _genomes or .est_)
    unwanted = [
        "genome",
        "genomic",
        "sequence",
        "dna",
        "cds",
        "pep",
        "transcript",
        "EST",
        "toplevel",
        "primary",
        "assembly",
    ]
    spacers = "( ?-?_?\.?)"  # noqa: W605
    for substring in unwanted:
        name = re.sub(
            f"{spacers}{substring}(s?){spacers}", "", name, flags=re.IGNORECASE
        )
    return name


def safe(name: Any) -> str:
    """Replace spaces with undescores."""
    return str(name).strip().replace(" ", "_")


def lower(string: Any) -> str:
    """safe(str).lower() for case-insensitive text comparisons"""
    return safe(string).lower()


def get_localname(name: Any, localname=None) -> str:
    """
    Returns the safe version of the given localname.

    If no localname is provided, return the safe version of the name.
    If the name is a working URL, return the safe version of the filename.
    """
    if localname:
        return safe(localname)

    # Remote file
    if try_except_pass((IOError, ValueError), urlopen, name):
        name = get_remotename(name)

    # Local file
    elif os.path.exists(name):
        name = get_genomename(name)

    return safe(name)


def try_except_pass(errors, func, *args, **kwargs):
    """
    try to return FUNC with ARGS, pass on ERRORS

    parameters
    ----------
    errors
      a single error, or a tuple of errors.

    func
      a function that takes args and kwargs
    """
    try:
        return func(*args, **kwargs)
    except errors:
        pass
