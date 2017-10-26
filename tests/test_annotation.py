from tempfile import mkdtemp, NamedTemporaryFile
import genomepy
import shutil
import gzip
import os

# Python 2
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

def setup():
    pass

def teardown():
    pass

def validate_gzipped_gtf(fname):
    assert os.path.exists(fname)
    with gzip.open(fname, 'r') as f:
        for line in f:
            line = line.decode()
            if line.startswith("#"):
                continue
            vals = line.split('\t')
            assert 9 == len(vals)
            int(vals[3]), int(vals[4])
            break

def validate_gzipped_bed(fname):
    assert os.path.exists(fname)
    with gzip.open(fname, 'r') as f:
        for line in f:
            line = line.decode()
            if line.startswith("#"):
                continue
            vals = line.split('\t')
            assert 12 == len(vals)
            int(vals[1]), int(vals[2])
            break

def test_annotation():
    tmp = mkdtemp()
    p = genomepy.provider.ProviderBase.create("UCSC")
    name = "sacCer3"
    
    p.download_annotation(name, tmp)
    
    gtf = os.path.join(tmp, name, name + ".annotation.gtf.gz")
    validate_gzipped_gtf(gtf)

    bed = os.path.join(tmp, name, name + ".annotation.bed.gz")
    validate_gzipped_bed(bed)
    
    shutil.rmtree(tmp)

def test_ensembl_annotation():
    tmp = mkdtemp()
    p = genomepy.provider.ProviderBase.create("Ensembl")
    
    for name, version in [("GRCh38.p10", 89), ("TAIR10", None)]:
        p.download_annotation(name, tmp)
    
        gtf = os.path.join(tmp, name, name + ".annotation.gtf.gz")
        validate_gzipped_gtf(gtf)
    
        bed = os.path.join(tmp, name, name + ".annotation.bed.gz")
        validate_gzipped_bed(bed)
    
    shutil.rmtree(tmp)


