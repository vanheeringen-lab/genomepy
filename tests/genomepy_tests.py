from nose.tools import *
from tempfile import mkdtemp
import genomepy
import shutil

def setup():
    print("SETUP!")

def teardown():
    print("TEAR DOWN!")

def test_basic():
    cfg = genomepy.functions.config
    assert 1 == len(cfg.keys())

@raises(FileNotFoundError)
def test_genome_dir_not_found():
    genomepy.genome("unknown", "unknown")

@raises(FileNotFoundError)
def test_no_fasta_files():
    genomepy.genome("empty", "tests/data/genome")

def test_ucsc_genome(): 
    """Test UCSC.
    
    Download S. cerevisiae genome from UCSC and retrieve a 
    specific sequence.
    """
    tmp = mkdtemp()
    genomepy.install_genome("sacCer3", "UCSC", genome_dir=tmp)
    g = genomepy.genome("sacCer3", genome_dir=tmp)
    seq = g["chrIV"][1337000:1337020] 
    assert str(seq) == "TTTGGTTGTTCCTCTTCCTT"
    shutil.rmtree(tmp)

def test_Ensembl_genome(): 
    """Test Ensembl.
    
    Download Drosophila genome from Ensembl and retrieve a 
    specific sequence.
    """
    tmp = mkdtemp()
    genomepy.install_genome("BDGP6", "Ensembl", genome_dir=tmp)
    g = genomepy.genome("BDGP6", genome_dir=tmp)
    seq = g["3L"][10637840:10637875] 
    assert str(seq).upper() == "TTTGCAACAGCTGCCGCAGTGTGACCGTTGTACTG"
    shutil.rmtree(tmp)
