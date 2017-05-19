from nose.tools import *
from tempfile import mkdtemp, NamedTemporaryFile
import genomepy
import shutil

# Python 2
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

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

def test_ensembl_genome(): 
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

def test_ncbi_genome(): 
    """Test NCBI.
    
    Download Drosophila genome from NCBI and retrieve a 
    specific sequence.
    """
    tmp = mkdtemp()
    genomepy.install_genome("Release 6 plus ISO1 MT", "NCBI", genome_dir=tmp)
    g = genomepy.genome("Release_6_plus_ISO1_MT", genome_dir=tmp)
    seq = g["3L"][10637840:10637875] 
    assert str(seq).upper() == "TTTGCAACAGCTGCCGCAGTGTGACCGTTGTACTG"
    shutil.rmtree(tmp)

def test_ucsc_human(): 
    """Test UCSC.
   
    Download human genome from UCSC and retrieve a 
    specific sequence.
    """
    tmp = mkdtemp()
    genomepy.install_genome("hg38", "UCSC", genome_dir=tmp)
    g = genomepy.genome("hg38", genome_dir=tmp)
    seq = g["chr6"][166168664:166168679] 
    assert str(seq) == "CCTCCTCGCTCTCTT"
    shutil.rmtree(tmp)

def test_ensembl_human(): 
    """Test Ensembl.
    
    Download human genome from Ensembl and retrieve a 
    specific sequence.
    """
    tmp = mkdtemp()
    genomepy.install_genome("GRCh38.p10", "Ensembl", genome_dir=tmp)
    g = genomepy.genome("GRCh38.p10", genome_dir=tmp)
    seq = g["6"][166168664:166168679] 
    assert str(seq) == "CCTCCTCGCTCTCTT"
    shutil.rmtree(tmp)

def test_ncbi_human(): 
    """Test NCBI.
    
    Download human genome from NCBI and retrieve a 
    specific sequence.
    """
    tmp = mkdtemp()
    genomepy.install_genome("GRCh38.p9", "NCBI", genome_dir=tmp)
    g = genomepy.genome("GRCh38.p9", genome_dir=tmp)
    seq = g["6"][166168664:166168679] 
    assert str(seq) == "CCTCCTCGCTCTCTT"
    shutil.rmtree(tmp)

def test_regexp_filter():
    fname = "tests/data/regexp/regexp.fa"

    regexps = [
        ('Chr.*', 2, 15),
        ('Scaffold.*', 1, 16),
        ('scaffold_.*', 3, 14),
        ('^\d+$', 4, 13),
        ('chr.*', 4, 13),
        ]

    tmpfa = NamedTemporaryFile(suffix=".fa").name
    for regex, match, no_match in regexps:
        fa = genomepy.utils.filter_fasta(
                fname, tmpfa, regex=regex, v=False, force=True)
        assert len(fa.keys()) == match
        fa = genomepy.utils.filter_fasta(
                fname, tmpfa, regex=regex, v=True, force=True)
        assert len(fa.keys()) == no_match

test_ncbi_human.slow = 1
test_ensembl_human.slow = 1
test_ucsc_human.slow = 1
