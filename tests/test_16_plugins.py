import os
import re
import subprocess as sp
from shutil import copyfile
from time import sleep

import pytest

import genomepy.utils
from genomepy.plugins import activate, deactivate, get_active_plugins, init_plugins
from genomepy.plugins.blacklist import BlacklistPlugin
from genomepy.plugins.bowtie2 import Bowtie2Plugin
from genomepy.plugins.bwa import BwaPlugin
from genomepy.plugins.gmap import GmapPlugin
from genomepy.plugins.hisat2 import Hisat2Plugin
from genomepy.plugins.minimap2 import Minimap2Plugin
from genomepy.plugins.star import StarPlugin
from tests import linux, travis


@pytest.fixture(autouse=True)
def activate_plugins():
    # save originally active plugins
    original_plugins = [p.name for p in get_active_plugins()]
    # activate all plugins
    [activate(p) for p in init_plugins()]

    yield

    # deactivate all plugins
    [deactivate(p) for p in init_plugins()]
    # reactivate original plugins
    [activate(p) for p in original_plugins]


@pytest.fixture(scope="function", params=["unzipped", "bgzipped"])
def genome(request):
    """Create a test genome and location"""
    name = "ce10"  # Use fake name for blacklist test
    fafile = "tests/data/small_genome.fa.gz"

    # setup test dir
    genomes_dir = os.path.join(os.getcwd(), ".genomepy_plugin_tests")
    genome_dir = os.path.join(genomes_dir, name)
    genomepy.utils.mkdir_p(genome_dir)
    fname = os.path.join(genome_dir, f"{name}.fa.gz")
    copyfile(fafile, fname)

    # unzip genome if required
    if request.param == "unzipped":
        sp.check_call(["gunzip", fname])

        # add annotation (for STAR and hisat2) for 1 of 2 tests
        gtf_file = "tests/data/ce10.annotation.gtf.gz"
        aname = os.path.join(genome_dir, f"{name}.annotation.gtf.gz")
        copyfile(gtf_file, aname)

    yield genomepy.Genome(name, genomes_dir=genomes_dir)

    # tear down test dir
    genomepy.utils.rm_rf(genomes_dir)


def dont_overwrite(p, genome, fname):
    t0 = os.path.getmtime(fname)
    # OSX rounds down getmtime to the second
    if not linux:
        sleep(1)
    p.after_genome_download(genome, force=False)
    t1 = os.path.getmtime(fname)
    assert t0 == t1


def test_blacklist(caplog, genome):
    """Create blacklist."""
    # independent of bgzipping
    # no need to check for both .fa and .fa.gz.
    if genome.filename.endswith(".fa.gz"):
        pass

    p = BlacklistPlugin()
    fname = re.sub(".fa(.gz)?$", ".blacklist.bed", genome.filename)

    # no blacklist found
    genome.name = "ce01"
    p.after_genome_download(genome, force=True)
    assert f"No blacklist found for {genome.name}" in caplog.text

    # error downloading blacklist
    genome.name = "this was a triumph"
    p.after_genome_download(genome, force=True)
    link = "I'm making a note here: 'Huge success'"
    assert f"Could not download blacklist file from {link}" in caplog.text

    # download blacklist
    genome.name = "ce10"
    p.after_genome_download(genome, force=True)
    assert "ce10-blacklist" in caplog.text
    assert os.path.exists(fname)
    with open(fname) as blacklist:
        for line in blacklist:
            assert line.startswith("chr")
            break
    # os.unlink(fname)

    # # download Ensembl/NCBI blacklist
    # genome.name = "GRCh38"
    # p.after_genome_download(genome, force=True)
    # assert "ENCFF356LFX/@@download/ENCFF356LFX.bed.gz" in caplog.text
    # with open(fname) as blacklist:
    #     for line in blacklist:
    #         assert not line.startswith("chr")
    #         break

    # don't overwrite
    dont_overwrite(p, genome, fname)
    os.unlink(fname)

    # genome.name = "ce10"


def test_bowtie2(genome, threads=2):
    """Create bowtie2 index."""
    # can work with bgzipped genomes natively,
    # no need to check for both .fa and .fa.gz.
    if genome.filename.endswith(".fa"):
        pass

    p = Bowtie2Plugin()
    p.after_genome_download(genome, threads=threads, force=True)

    dirname = os.path.dirname(genome.filename)
    index_dir = os.path.join(dirname, "index", "bowtie2")
    fname = os.path.join(index_dir, f"{genome.name}.1.bt2")
    assert os.path.exists(index_dir)
    assert os.path.exists(fname)


def test_bwa(genome, threads=2):
    """Create bwa index."""
    # can work with bgzipped genomes natively,
    # no need to check for both .fa and .fa.gz.
    if genome.filename.endswith(".fa"):
        pass

    p = BwaPlugin()
    p.after_genome_download(genome, threads=threads, force=True)

    dirname = os.path.dirname(genome.filename)
    index_dir = os.path.join(dirname, "index", "bwa")
    fname = os.path.join(index_dir, f"{genome.name}.fa.sa")
    assert os.path.exists(index_dir)
    assert os.path.exists(fname)


@pytest.mark.skipif(not travis, reason="slow")
def test_gmap(genome, threads=2):
    """Create gmap index."""
    p = GmapPlugin()
    p.after_genome_download(genome, threads=threads, force=True)

    dirname = os.path.dirname(genome.filename)
    index_dir = os.path.join(dirname, "index", "gmap")
    fname = os.path.join(index_dir, f"{genome.name}.maps")
    assert os.path.exists(index_dir)
    assert os.path.exists(fname)


def test_hisat2(caplog, genome, threads=2):
    """Create hisat2 index."""
    p = Hisat2Plugin()
    p.after_genome_download(genome, threads=threads, force=True)

    dirname = os.path.dirname(genome.filename)
    index_dir = os.path.join(dirname, "index", "hisat2")
    fname = os.path.join(index_dir, f"{genome.name}.1.ht2")
    assert os.path.exists(index_dir)
    assert os.path.exists(fname)

    if genome.annotation_gtf_file:
        # check if splice-aware index is generated
        assert os.path.exists(os.path.join(genome.genome_dir, "splice_sites.txt"))
        assert os.path.exists(os.path.join(genome.genome_dir, "exon_sites.txt"))
        # check if annotation file is still the same
        assert os.path.exists(genome.annotation_gtf_file)
        assert genome.annotation_gtf_file.endswith(".gtf.gz")
    else:
        assert "Creating Hisat2 index without annotation file." in caplog.text


def test_minimap2(genome, threads=2):
    """Create minimap2 index."""
    # can work with bgzipped genomes natively,
    # no need to check for both .fa and .fa.gz.
    if genome.filename.endswith(".fa"):
        pass

    p = Minimap2Plugin()
    p.after_genome_download(genome, threads=threads, force=True)

    dirname = os.path.dirname(genome.filename)
    index_dir = os.path.join(dirname, "index", "minimap2")
    fname = os.path.join(index_dir, f"{genome.name}.mmi")
    assert os.path.exists(index_dir)
    assert os.path.exists(fname)


@pytest.mark.skipif(not travis, reason="slow")
def test_star(caplog, genome, threads=2):
    """Create star index."""
    p = StarPlugin()
    p.after_genome_download(genome, threads=threads, force=True)

    dirname = os.path.dirname(genome.filename)
    index_dir = os.path.join(dirname, "index", "star")
    fname = os.path.join(index_dir, "SA")
    assert os.path.exists(index_dir)
    assert os.path.exists(fname)

    if genome.annotation_gtf_file:
        # check if splice-aware index is generated
        assert "Creating star index..." in caplog.text
        # check if annotation file is still the same
        assert os.path.exists(genome.annotation_gtf_file)
        assert genome.annotation_gtf_file.endswith(".gtf.gz")
    else:
        assert "Creating STAR index without annotation file." in caplog.text
