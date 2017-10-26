from tempfile import mkdtemp, NamedTemporaryFile
from genomepy.utils import generate_gap_bed
import shutil
import gzip
import os

# Python 2
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

def test_gaps():
    infa = "tests/data/gap.fa"
    outbed = "tests/data/gap.bed" 

    tmp = NamedTemporaryFile().name
    generate_gap_bed(infa, tmp)

    result = open(tmp).read()
    expect = open(outbed).read()
    
    assert result == expect
