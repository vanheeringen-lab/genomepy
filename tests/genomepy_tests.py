from nose.tools import *
import genomepy

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
