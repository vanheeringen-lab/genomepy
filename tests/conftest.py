"""
Global fixtures and functions for pytest
pytest can only share fixtures between modules if they are declared here.
"""
import logging
import os

import pytest
from loguru import logger

import genomepy.providers
from genomepy.providers.base import BaseProvider
from genomepy.providers.ensembl import EnsemblProvider
from genomepy.providers.gencode import GencodeProvider
from genomepy.providers.local import LocalProvider
from genomepy.providers.ncbi import NcbiProvider
from genomepy.providers.ucsc import UcscProvider
from genomepy.providers.url import UrlProvider


@pytest.fixture(scope="function")
def caplog(caplog):
    """Fixture is necessary to be able to check loguru log messages"""

    class PropogateHandler(logging.Handler):
        def emit(self, record):
            logging.getLogger(record.name).handle(record)

    handler_id = logger.add(PropogateHandler(), format="{message} {extra}")
    yield caplog
    logger.remove(handler_id)


def teardown(gprefix, skip=None):
    for ext in [
        "fa.fai",
        "fa.sizes",
        "gaps.bed",
        "fa.gz.fai",
        "fa.gz.sizes",
        "annotation.gtf",
        "annotation.bed",
    ]:
        if skip and ext in skip:
            continue
        file = gprefix + ext
        if os.path.exists(file):
            os.remove(file)
    gdir = os.path.dirname(gprefix)
    readme = os.path.join(gdir, "README.txt")
    if os.path.exists(readme):
        os.remove(readme)


@pytest.fixture(scope="function")
def small_genome():
    yield genomepy.Genome("tests/data/small_genome.fa.gz")

    teardown("tests/data/small_genome.")


@pytest.fixture(scope="function")
def gap_genome():
    yield genomepy.Genome("tests/data/gap.fa")

    teardown("tests/data/gap.")


@pytest.fixture(scope="function")
def annot():
    genome_file = "tests/data/regexp/regexp.fa"
    gtf_file = "tests/data/regexp/regexp.annotation.gtf"
    bed_file = "tests/data/regexp/regexp.annotation.bed"

    genomepy.Genome(genome_file)
    with open(gtf_file, "w") as f:
        f.write("# skip this line\n")
        f.write(
            """chrM\tvanHeeringen-lab\tNP_059343.1\t15307\t16448\t42\t+\t.\tattributes"""
        )
    with open(bed_file, "w") as f:
        f.write(
            """chrM\t15307\t16448\tNP_059343.1\t42\t+\t15307\t16448\t0\t1\t1141,\t0,"""
        )
    yield genomepy.Annotation("regexp", genomes_dir="tests/data")

    teardown("tests/data/regexp/regexp.")


def validate_annot(fname, ftype):
    """fname = path, ftype = 'bed' or 'gtf'."""
    assert os.path.exists(fname)
    columns = 12 if ftype == "bed" else 9
    start, end = (3, 4) if ftype == "gtf" else (1, 2)
    with open(fname, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            vals = line.split("\t")
            assert columns == len(vals)
            int(vals[start]), int(vals[end])
            break


@pytest.fixture(scope="function")
def base():
    return BaseProvider()


@pytest.fixture(scope="function")
def ensembl():
    return EnsemblProvider()


@pytest.fixture(scope="function")
def ucsc():
    return UcscProvider()


@pytest.fixture(scope="function")
def gencode():
    return GencodeProvider()


@pytest.fixture(scope="function")
def ncbi():
    return NcbiProvider()


@pytest.fixture(scope="function")
def local():
    return LocalProvider()


@pytest.fixture(scope="function")
def url():
    return UrlProvider()


@pytest.fixture(scope="function")
def provider():
    return genomepy.Provider()
