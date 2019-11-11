import os
import pytest
import re
from subprocess import check_call
from tempfile import mkdtemp
from shutil import rmtree, copyfile

from genomepy.plugin import init_plugins, activate
from genomepy.utils import cmd_ok
from genomepy.utils import mkdir_p
from genomepy.functions import Genome
from genomepy.plugins.bwa import BwaPlugin
# from genomepy.plugins.gmap import GmapPlugin
from genomepy.plugins.minimap2 import Minimap2Plugin
from genomepy.plugins.bowtie2 import Bowtie2Plugin
from genomepy.plugins.hisat2 import Hisat2Plugin
from genomepy.plugins.blacklist import BlacklistPlugin


@pytest.fixture(scope="module")
def tempdir():
    """Temporary directory."""
    tmpdir = mkdtemp()
    yield tmpdir
    rmtree(tmpdir)


@pytest.fixture(scope="module", params=["unzipped", "bgzipped"])
def genome(request, tempdir):
    """Create a test genome."""
    name = "dm3"  # Use fake name for blacklist test
    fafile = "tests/data/small_genome.fa"
    bgzipped = True if request.param == "bgzipped" else False

    # Input needs to be bgzipped, depending on param
    if os.path.exists(fafile + ".gz"):
        if not bgzipped:
            check_call(["gunzip", fafile + ".gz"])
    elif bgzipped:
        check_call(["bgzip", fafile])

    tmpdir = os.path.join(tempdir, request.param, name)
    mkdir_p(tmpdir)

    if bgzipped:
        fafile += ".gz"

    copyfile(fafile, os.path.join(tmpdir, os.path.basename(fafile)))
    for p in init_plugins():
        activate(p)
    # provide the fixture value
    yield Genome(name, genome_dir=os.path.join(tempdir, request.param))
    if os.path.exists(fafile) and not bgzipped:
        check_call(["bgzip", fafile])


@pytest.fixture(scope="module", params=["no-overwrite", "overwrite"])
def force(request):
    return request.param


def test_blacklist(genome, force):
    """Create blacklist."""
    assert os.path.exists(genome.filename)

    force = True if force == "overwrite" else False
    p = BlacklistPlugin()

    p.after_genome_download(genome, force=force)
    fname = re.sub(".fa(.gz)?$", ".blacklist.bed.gz", genome.filename)
    assert os.path.exists(fname)

    t0 = os.path.getmtime(fname)
    import sys
    sys.stderr.write(str(t0))
    p.after_genome_download(genome, force=force)
    t1 = os.path.getmtime(fname)
    sys.stderr.write(str(t1))
    sys.stderr.write(str(type(t1)))
    assert t0 != t1 if force else t0 == t1


def test_bwa(genome, force):
    """Create bwa index."""
    assert os.path.exists(genome.filename)

    force = True if force == "overwrite" else False
    if cmd_ok("bwa"):
        p = BwaPlugin()
        p.after_genome_download(genome, force=force)
        dirname = os.path.dirname(genome.filename)
        index_dir = os.path.join(dirname, "index", "bwa")
        fname = os.path.join(index_dir, "{}.fa.sa".format(genome.name))
        assert os.path.exists(index_dir)
        assert os.path.exists(fname)

        t0 = os.path.getmtime(fname)
        p.after_genome_download(genome, force=force)
        t1 = os.path.getmtime(fname)
        assert t0 != t1 if force else t0 == t1


def test_minimap2(genome, force):
    """Create minimap2 index."""
    assert os.path.exists(genome.filename)

    force = True if force == "overwrite" else False
    if cmd_ok("minimap2"):
        p = Minimap2Plugin()
        p.after_genome_download(genome, force=force)
        dirname = os.path.dirname(genome.filename)
        index_dir = os.path.join(dirname, "index", "minimap2")
        fname = os.path.join(index_dir, "{}.mmi".format(genome.name))
        assert os.path.exists(index_dir)
        assert os.path.exists(fname)

        t0 = os.path.getmtime(fname)
        p.after_genome_download(genome, force=force)
        t1 = os.path.getmtime(fname)
        assert t0 != t1 if force else t0 == t1


def test_bowtie2(genome, force):
    """Create bbowtie2 index."""
    assert os.path.exists(genome.filename)

    force = True if force == "overwrite" else False
    if cmd_ok("bowtie2"):
        p = Bowtie2Plugin()
        p.after_genome_download(genome, force=force)
        dirname = os.path.dirname(genome.filename)
        index_dir = os.path.join(dirname, "index", "bowtie2")
        fname = os.path.join(index_dir, "{}.1.bt2".format(genome.name))
        assert os.path.exists(index_dir)
        assert os.path.exists(fname)

        t0 = os.path.getmtime(fname)
        p.after_genome_download(genome, force=force)
        t1 = os.path.getmtime(fname)
        assert t0 != t1 if force else t0 == t1


def test_hisat2(genome, force):
    """Create hisat2 index."""
    assert os.path.exists(genome.filename)

    force = True if force == "overwrite" else False
    if cmd_ok("hisat2-build"):
        p = Hisat2Plugin()
        p.after_genome_download(genome)
        dirname = os.path.dirname(genome.filename)
        index_dir = os.path.join(dirname, "index", "hisat2")
        fname = os.path.join(index_dir, "{}.1.ht2".format(genome.name))
        assert os.path.exists(index_dir)
        assert os.path.exists(fname)

        t0 = os.path.getmtime(fname)
        p.after_genome_download(genome, force=force)
        t1 = os.path.getmtime(fname)
        assert t0 != t1 if force else t0 == t1


# def test_gmap(genome, force):
#     """Create gmap index."""
#     assert os.path.exists(genome.filename)
#
#     force = True if force == "overwrite" else False
#     if cmd_ok("gmap"):
#         p = GmapPlugin()
#         p.after_genome_download(genome, force=force)
#         dirname = os.path.dirname(genome.filename)
#         index_dir = os.path.join(dirname, "index", "gmap")
#         fname = os.path.join(index_dir, "{}.maps".format(genome.name))
#         assert os.path.exists(index_dir)
#         assert os.path.exists(fname)
#
#         t0 = os.path.getmtime(fname)
#         p.after_genome_download(genome, force=force)
#         t1 = os.path.getmtime(fname)
#         assert t0 != t1 if force else t0 == t1
