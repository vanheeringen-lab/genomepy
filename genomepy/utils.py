"""Utility functions."""
import os
import itertools
import re
import sys
import subprocess as sp
import shutil
import time
from typing import Optional, List
from urllib.request import urlopen

from loguru import logger
import norns
from norns import exceptions

config = norns.config("genomepy", default="cfg/default.yaml")


def mkdir_p(path):
    """'mkdir -p' in Python"""
    path = os.path.expanduser(path)
    os.makedirs(path, exist_ok=True)


def rm_rf(path):
    """'rm -rf' in Python"""
    path = os.path.expanduser(path)
    if os.path.isfile(path):
        try:
            os.unlink(path)
        except OSError:  # in case of NTFS related issues
            pass
    elif os.path.isdir(path):
        shutil.rmtree(path, ignore_errors=True)


def cmd_ok(cmd):
    """Returns True if cmd can be run."""
    try:
        sp.check_call(cmd, stderr=sp.PIPE, stdout=sp.PIPE)
    except sp.CalledProcessError:
        # bwa gives return code of 1 with no argument
        pass
    except FileNotFoundError:
        logger.error(f"{cmd} not found, skipping\n")
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


def get_genomes_dir(genomes_dir: str = None, check_exist: Optional[bool] = True) -> str:
    """import genomes_dir if none is given, and check validity"""
    if not genomes_dir:
        # backwards compatibility for "genome_dir" (this fixes issue #87)
        genomes_dir = config.get("genomes_dir", config.get("genome_dir", None))
    if not genomes_dir:
        raise exceptions.ConfigError("Please provide or configure a genomes_dir")

    genomes_dir = os.path.abspath(os.path.expanduser(genomes_dir))
    if not os.path.exists(genomes_dir) and check_exist:
        raise FileNotFoundError(f"Genomes_dir {genomes_dir} does not exist!")

    return genomes_dir


def safe(name) -> str:
    """Replace spaces with undescores."""
    return str(name).strip().replace(" ", "_")


def lower(string) -> str:
    """for case-insensitive text comparisons"""
    return safe(string).lower()


def get_localname(name, localname=None):
    """
    Returns the safe version of the given localname, if provided.
    If no localname is provided, return the safe version of the name.
    If the name is a working URL, return the safe version of the filename.
    """
    if localname:
        return safe(localname)
    try:
        urlopen(name)
    except (IOError, ValueError):
        return safe(name)
    else:
        # try to get the name from the url
        name = name.split("/")[-1]  # remove path
        name = name.replace(".gz", "")  # remove .gz
        name = os.path.splitext(name)[0]  # remove .fa/.fna/.fasta etc
        name = safe(name)  # remove spaces
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
        spacers = "(-?_?\.?)"  # noqa: W605
        for substring in unwanted:
            name = re.sub(
                f"{spacers}{substring}(s?){spacers}", "", name, flags=re.IGNORECASE
            )
        return name


def try_except_pass(errors, func, *args, **kwargs):
    """try to return FUNC with ARGS, pass on ERRORS"""
    try:
        return func(*args, **kwargs)
    except errors:
        pass


def best_accession(reference: str, targets: list):
    """Return the nearest accession ID from a list of IDs"""
    if len(targets) == 1:
        return targets[0]

    # GCA/GCF
    matching_prefix = [t for t in targets if t.startswith(reference[0:3])]
    if len(matching_prefix) > 0:
        targets = matching_prefix

    # patch levels
    # e.g. GCA_000002035.4 & GCA_000002035.3
    ref_patch = int(reference.split(".")[1]) if "." in reference else 0
    nearest = [999, []]
    for target in targets:
        tgt_patch = int(target.split(".")[1]) if "." in target else 0
        distance = abs(ref_patch - tgt_patch)
        if distance == nearest[0]:
            nearest[1].append(target)
        if distance < nearest[0]:
            nearest = [distance, [target]]
    targets = nearest[1]

    if len(targets) > 1:
        logger.warning(
            f"Multiple matching accession numbers found. Returning the closest ({targets[0]})."
        )
    return targets[0]


def best_search_result(asm_acc: str, results: List[list]) -> list:
    """Return the best result from ProviderBase.search_all based on accession IDs"""
    if len(results) == 0:
        logger.warning(f"No assembly found similar to {asm_acc}")
        return []

    if len(results) > 1:
        accessions = [res[2] for res in results]
        bes_acc = best_accession(asm_acc, accessions)
        results = [r for r in results if r[2] == bes_acc]

    return results[0]
