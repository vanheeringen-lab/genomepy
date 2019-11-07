from tempfile import mkdtemp, NamedTemporaryFile
import genomepy
import shutil
import pytest
import os

# Python 2
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"


def test_basic():
    cfg = genomepy.functions.config
    print(cfg)
    assert 3 == len(cfg.keys())


def test_genome_dir_not_found():
    with pytest.raises(FileNotFoundError):
        genomepy.Genome("unknown", "unknown")


def test_no_fasta_files():
    with pytest.raises(FileNotFoundError):
        genomepy.Genome("empty", "tests/data/genome")


def test_track_type():
    tracks = [
        (("chr1:10-20", "chr2:10-20"), "interval"),
        (["chr1:10-20", "chr2:10-20"], "interval"),
        ("tests/data/regions.txt", "interval"),
        ("tests/data/regions.bed", "bed"),
    ]

    for track, track_type in tracks:
        result = genomepy.functions.get_track_type(track)
        assert result == track_type


@pytest.fixture(
    scope="module",
    params=[
        ["unzipped", "no-overwrite"],
        ["unzipped", "overwrite"],
        ["bgzipped", "no-overwrite"],
        ["bgzipped", "overwrite"],
    ],
)
def status(request):
    return [request.param[0], request.param[1]]


def test_install_genome(status, genome="sacCer3"):
    """Test UCSC, force and bgzip

    Download S. cerevisiae genome from UCSC and retrieve a specific sequence.

    Force and bgzip are called regardless of provider, and are therefore only tested here
    """
    tmp = mkdtemp()
    bgzip = False
    if status[0] == "bgzipped":
        bgzip = True
    force = False
    if status[1] == "overwrite":
        force = True

    genomepy.install_genome(genome, "UCSC", genome_dir=tmp, bgzip=bgzip, force=force)
    g = genomepy.Genome("sacCer3", genome_dir=tmp)
    seq = g["chrIV"][1337000:1337020]
    assert str(seq) == "TTTGGTTGTTCCTCTTCCTT"

    # Check if force==True overwrites files, and if force==False doesn't
    ext = ".fa.gz" if bgzip else ".fa"
    path = os.path.join(tmp, genome, genome + ext)
    t0 = os.path.getmtime(path)
    genomepy.install_genome(genome, "UCSC", genome_dir=tmp, bgzip=bgzip, force=force)
    t1 = os.path.getmtime(path)
    if force:
        assert t0 != t1
    else:
        assert t0 == t1

    shutil.rmtree(tmp)


def test_ucsc_genome():
    """Already tested ucsc in previous test`"""
    assert 1


# 2019-11-07 BDGP6, BDGP6.22 and dere_caf1 currently fails on Ensembl
@pytest.mark.xfail()
def test_ensembl_genome(genome="dere_caf1"):
    """Test Ensembl.

    Download Drosophila genome from Ensembl and retrieve a
    specific sequence.
    """
    tmp = mkdtemp()
    genomepy.install_genome(genome, "Ensembl", genome_dir=tmp)
    g = genomepy.Genome(genome, genome_dir=tmp)
    seq = g["scaffold_4929"][61:121]
    assert (
        str(seq).upper()
        == "ACACTTGGCTAAATTCCAAGAGAACTTACAAAAACAAACACCCCCCAGTCGAAAAAGGCC"
    )
    # seq = g["3L"][10637840:10637875]
    # assert str(seq).upper() == "TTTGCAACAGCTGCCGCAGTGTGACCGTTGTACTG"
    shutil.rmtree(tmp)


def test_ensembl_genome_download_links():
    """Test Ensembl FTP links for various genomes

    These genomes are hosted on ftp.ensembl.org.
    """
    p = genomepy.provider.ProviderBase.create("Ensembl")

    for genome in ["GRCz11", "GRCh38.p13"]:
        p.get_genome_download_link(genome)


@pytest.mark.skipif(travis, reason="FTP does not work on Travis")
def test_ensemblgenomes_genome_download_links():
    """Test Ensembl FTP links for various genomes

    These genomes are hosted on ftp.ensemblgenomes.org.
    """
    p = genomepy.provider.ProviderBase.create("Ensembl")

    for genome in ["AgamP4", "WBcel235", "R64-1-1"]:
        p.get_genome_download_link(genome)


def test_ncbi_genome():
    """Test NCBI.

    Download Drosophila genome from NCBI and retrieve a
    specific sequence.
    """
    tmp = mkdtemp()
    genomepy.install_genome("Release 6 plus ISO1 MT", "NCBI", genome_dir=tmp)
    g = genomepy.Genome("Release_6_plus_ISO1_MT", genome_dir=tmp)
    seq = g["3L"][10637840:10637875]
    assert str(seq).upper() == "TTTGCAACAGCTGCCGCAGTGTGACCGTTGTACTG"
    shutil.rmtree(tmp)


def test_regexp_filter():
    fname = "tests/data/regexp/regexp.fa"

    regexps = [
        ("Chr.*", 2, 15),
        ("Scaffold.*", 1, 16),
        ("scaffold_.*", 3, 14),
        (r"^\d+$", 4, 13),
        ("chr.*", 4, 13),
    ]

    tmpfa = NamedTemporaryFile(suffix=".fa").name
    for regex, match, no_match in regexps:
        fa = genomepy.utils.filter_fasta(fname, tmpfa, regex=regex, v=False, force=True)
        assert len(fa.keys()) == match
        fa = genomepy.utils.filter_fasta(fname, tmpfa, regex=regex, v=True, force=True)
        assert len(fa.keys()) == no_match


def test_url_genome():
    """Test URL.

    Download S. cerevisiae genome directly from an url from UCSC and retrieve a
    specific sequence.
    """
    tmp = mkdtemp()
    genomepy.install_genome(
        "http://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/chromFa.tar.gz",
        "url",
        genome_dir=tmp,
        localname="url_test",
    )
    g = genomepy.Genome("url_test", genome_dir=tmp)
    assert str(g["chrI"][:12]).lower() == "gcctaagcctaa"
    shutil.rmtree(tmp)


def test_bad_url():
    """Test URL.

    Try to download a non-genome url link.
    """
    with pytest.raises(ValueError):
        genomepy.install_genome("www.google.com", "url")


def test_nonexisting_url():
    """Test URL.

    Try to download a non-genome url link.
    """
    with pytest.raises(ValueError):
        genomepy.install_genome("this is not an url", "url")
