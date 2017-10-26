import os
import pytest
from subprocess import check_call
from tempfile import mkdtemp
from shutil import rmtree, copyfile

from genomepy.plugin import init_plugins, activate
from genomepy.functions import Genome
from genomepy.plugins.bwa import BwaPlugin
from genomepy.plugins.gmap import GmapPlugin
from genomepy.plugins.minimap2 import Minimap2Plugin
from genomepy.plugins.bowtie2 import Bowtie2Plugin
from genomepy.plugins.hisat2 import Hisat2Plugin

@pytest.fixture(scope="module")
def tempdir():
    """Temporary directory.""" 
    tmpdir = mkdtemp()
    yield tmpdir
    rmtree(tmpdir)

@pytest.fixture(scope="module")
def genome(tempdir):
    """Create a test genome.""" 
    name = "small_genome" 
    fafile = "tests/data/small_genome.fa"
    if os.path.exists(fafile + ".gz"):
        check_call(["gunzip", fafile + ".gz"])

    os.mkdir(os.path.join(tempdir, name))
    copyfile(fafile, os.path.join(tempdir, name, os.path.basename(fafile)))
    for p in init_plugins():
        activate(p)
    yield Genome(name, genome_dir=tempdir)  # provide the fixture value
    if os.path.exists(fafile):
        check_call(["gzip", fafile])

def test_bwa(genome):
    """Create bwa index.""" 
    assert os.path.exists(genome.filename)
    p = BwaPlugin()
    p.after_genome_download(genome)
    dirname = os.path.dirname(genome.filename)
    index_dir = os.path.join(dirname, "index" , "bwa")
    assert os.path.exists(index_dir)
    assert os.path.exists(os.path.join(index_dir, "{}.fa.sa".format(genome.name)))

def test_minimap2(genome):
    """Create minimap2 index.""" 
    assert os.path.exists(genome.filename)
    p = Minimap2Plugin()
    p.after_genome_download(genome)
    dirname = os.path.dirname(genome.filename)
    index_dir = os.path.join(dirname, "index" , "minimap2")
    assert os.path.exists(index_dir)
    assert os.path.exists(os.path.join(index_dir, "{}.mmi".format(genome.name)))

def test_bowtie2(genome):
    """Create bowtie2 index.""" 
    assert os.path.exists(genome.filename)
    p = Bowtie2Plugin()
    p.after_genome_download(genome)
    dirname = os.path.dirname(genome.filename)
    index_dir = os.path.join(dirname, "index" , "bowtie2")
    assert os.path.exists(index_dir)
    assert os.path.exists(os.path.join(index_dir, "{}.1.bt2".format(genome.name)))

def test_hisat2(genome):
    """Create hisat2 index.""" 
    assert os.path.exists(genome.filename)
    p = Hisat2Plugin()
    p.after_genome_download(genome)
    dirname = os.path.dirname(genome.filename)
    index_dir = os.path.join(dirname, "index" , "hisat2")
    assert os.path.exists(index_dir)
    assert os.path.exists(os.path.join(index_dir, "{}.1.ht2".format(genome.name)))

def test_gmap(genome):
    """Create gmap index.""" 
    assert os.path.exists(genome.filename)
    p = GmapPlugin()
    p.after_genome_download(genome)
    dirname = os.path.dirname(genome.filename)
    index_dir = os.path.join(dirname, "index" , "gmap", genome.name)
    assert os.path.exists(index_dir)
    assert os.path.exists(os.path.join(index_dir, "{}.version".format(genome.name)))
