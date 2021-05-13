"""
Global fixtures and functions for pytest
pytest can only share fixtures between modules if they are declared here.
"""
import os
import logging

from loguru import logger
import pytest

import genomepy


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
    yield genomepy.Annotation(name="regexp", genomes_dir="tests/data")

    teardown("tests/data/regexp/regexp.", skip=["fa.fai"])


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


@pytest.fixture(scope="package")
def provider():
    return genomepy.ProviderBase()


@pytest.fixture(scope="package")
def ensembl(provider):
    return provider.create("ensembl")


@pytest.fixture(scope="package")
def ucsc(provider):
    return provider.create("ucsc")


@pytest.fixture(scope="package")
def ncbi(provider):
    return provider.create("ncbi")


@pytest.fixture(scope="package")
def url(provider):
    return provider.create("url")


# # turns out to be slower than linear processing...
#
# from multiprocessing.pool import Pool
#
# # start downloading provider genomes in the background
# pool = Pool(processes=1)
# async_ncbi = pool.apply_async(func=genomepy.ProviderBase.create, kwds={"name": "ncbi"})
# async_ucsc = pool.apply_async(func=genomepy.ProviderBase.create, kwds={"name": "ucsc"})
# async_ens = pool.apply_async(
#     func=genomepy.ProviderBase.create, kwds={"name": "ensembl"}
# )
# pool.close()
#
#
# @pytest.fixture(scope="package")
# def provider():
#     return genomepy.provider.ProviderBase()
#
#
# @pytest.fixture(scope="package")
# def ensembl():
#     return async_ens.get()
#
#
# @pytest.fixture(scope="package")
# def ucsc():
#     return async_ucsc.get()
#
#
# @pytest.fixture(scope="package")
# def ncbi():
#     return async_ncbi.get()
#
#
# @pytest.fixture(scope="package")
# def url(provider):
#     return provider.create("url")
