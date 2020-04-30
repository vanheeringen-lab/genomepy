import os
import genomepy
import pytest
import re
import shutil
import subprocess as sp
from platform import system

from genomepy.plugin import init_plugins, activate, deactivate
from genomepy.plugins.blacklist import BlacklistPlugin
from genomepy.plugins.bowtie2 import Bowtie2Plugin
from genomepy.plugins.bwa import BwaPlugin
from genomepy.plugins.gmap import GmapPlugin
from genomepy.plugins.hisat2 import Hisat2Plugin
from genomepy.plugins.minimap2 import Minimap2Plugin
from genomepy.plugins.star import StarPlugin

linux = system() == "Linux"
travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"


def test_plugins():
    # activate and check all plugins
    for p in init_plugins():
        if p not in ["blacklist", "star"]:
            assert genomepy.utils.cmd_ok(p)
        elif p == "star":
            assert genomepy.utils.cmd_ok(p.upper())
        activate(p)


@pytest.fixture(scope="module", params=["unzipped", "bgzipped"])
def genome(request):
    """Create a test genome and location"""
    name = "ce10"  # Use fake name for blacklist test
    fafile = "tests/data/small_genome.fa.gz"

    genomes_dir = os.path.join(os.getcwd(), ".genomepy_plugin_tests")
    if os.path.exists(genomes_dir):
        shutil.rmtree(genomes_dir)
    genome_dir = os.path.join(genomes_dir, name)
    genomepy.utils.mkdir_p(genome_dir)
    fname = os.path.join(genome_dir, f"{name}.fa.gz")
    shutil.copyfile(fafile, fname)

    # unzip genome if required
    if request.param == "unzipped":
        sp.check_call(["gunzip", fname])

        # add annotation (for STAR and hisat2), but only once
        gtf_file = "tests/data/ce10.annotation.gtf"
        aname = os.path.join(genome_dir, f"{name}.annotation.gtf")
        shutil.copyfile(gtf_file, aname)

    return genomepy.Genome(name, genomes_dir=genomes_dir)


def test_blacklist(genome):
    """Create blacklist."""
    # not affected by bgzipping,
    # no need to check for both .fa and .fa.gz.
    if genome.filename.endswith(".fa"):
        pass

    p = BlacklistPlugin()
    p.after_genome_download(genome)

    fname = re.sub(".fa(.gz)?$", ".blacklist.bed", genome.filename)
    assert os.path.exists(fname)


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


@pytest.mark.skipif(not travis or not linux, reason="slow")
def test_gmap(genome, threads=2):
    """Create gmap index."""
    p = GmapPlugin()
    p.after_genome_download(genome, threads=threads, force=True)

    dirname = os.path.dirname(genome.filename)
    index_dir = os.path.join(dirname, "index", "gmap")
    fname = os.path.join(index_dir, f"{genome.name}.maps")
    assert os.path.exists(index_dir)
    assert os.path.exists(fname)


def test_hisat2(capsys, genome, threads=2):
    """Create hisat2 index."""
    p = Hisat2Plugin()
    p.after_genome_download(genome, threads=threads, force=True)

    dirname = os.path.dirname(genome.filename)
    index_dir = os.path.join(dirname, "index", "hisat2")
    fname = os.path.join(index_dir, f"{genome.name}.1.ht2")
    assert os.path.exists(index_dir)
    assert os.path.exists(fname)

    # check if splice-aware index is generated
    captured = capsys.readouterr().out.strip()
    if genome.annotation_gtf_file:
        assert captured == ""
    else:
        assert captured == "Generating Hisat2 index without annotation file."


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


@pytest.mark.skipif(not travis or not linux, reason="slow")
def test_star(capsys, genome, threads=2):
    """Create star index."""
    p = StarPlugin()
    p.after_genome_download(genome, threads=threads, force=True)

    dirname = os.path.dirname(genome.filename)
    index_dir = os.path.join(dirname, "index", "star")
    fname = os.path.join(index_dir, "SA")
    assert os.path.exists(index_dir)
    assert os.path.exists(fname)

    # check if splice-aware index is generated
    captured = capsys.readouterr().out.strip()
    if genome.annotation_gtf_file:
        assert captured == ""
    else:
        assert captured == "Generating STAR index without annotation file."


def test_plugin_cleanup():
    for p in init_plugins():
        deactivate(p)

    # cleanup after testing pluging
    genome_dir = os.path.join(os.getcwd(), ".genomepy_plugin_tests")
    shutil.rmtree(genome_dir)
