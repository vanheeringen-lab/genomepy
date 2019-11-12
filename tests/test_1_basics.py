import genomepy
import pytest
import os

# # Python 2
# try:
#     FileNotFoundError
# except NameError:
#     FileNotFoundError = IOError

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
