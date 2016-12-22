from nose.tools import *
import genomepy

def setup():
    print "SETUP!"

def teardown():
    print "TEAR DOWN!"

def test_basic():
    cfg = genomepy.config
    assert 1 == len(cfg.keys())
