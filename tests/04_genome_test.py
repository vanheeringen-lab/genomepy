import genomepy
import os
import pytest


# to ignore file changes
# git update-index --assume-unchanged tests/data/small_genome.fa.gz
# to recheck file changes
# git update-index --no-assume-unchanged tests/data/small_genome.fa.gz
def test_genome__init__(genome="tests/data/small_genome.fa.gz"):
    # no fasta file
    with pytest.raises(FileNotFoundError):
        genomepy.Genome("empty", "tests/data/genome")

    # genome dir not found
    with pytest.raises(FileNotFoundError):
        genomepy.Genome("unknown", "unknown")

    # initialize the class (creates the index file)
    index_file = genome + ".fai"
    if os.path.exists(index_file):
        os.unlink(index_file)

    genomepy.Genome(genome)
    assert os.path.exists(index_file)
    os.unlink(index_file)  # TODO: remove?


def test__read_metadata(genome="tests/data/small_genome.fa.gz"):
    pass


def test__bed_to_seqs(
    genome="tests/data/small_genome.fa.gz", track="tests/data/regions.bed"
):
    g = genomepy.Genome(genome)

    # extract sequences marked in regions.bed from small_genome.fa.gz
    seqs = g._bed_to_seqs(track=track, stranded=False, extend_up=0, extend_down=0)
    for i, seq in enumerate(seqs):
        assert seq.name == ["chrI:10-20 gene_a", "chrII:20-30 gene_b"][i]
        assert seq.seq == ["CCCACACACC", "TCCTCCAAGC"][i]

    # second sequence is on the negative strand
    seqs = g._bed_to_seqs(track=track, stranded=True, extend_up=0, extend_down=0)
    for i, seq in enumerate(seqs):
        assert seq.name == ["chrI:10-20 gene_a", "chrII:20-30 gene_b"][i]
        # original:        "CCCACACACC", "TCCTCCAAGC"
        assert seq.seq == ["CCCACACACC", "GCTTGGAGGA"][i]

    # extend by varying amounts
    seqs = g._bed_to_seqs(track=track, stranded=True, extend_up=1, extend_down=2)
    for i, seq in enumerate(seqs):
        assert seq.name == ["chrI:10-20 gene_a", "chrII:20-30 gene_b"][i]
        # original:         "CCCACACACC",    "GCTTGGAGGA"
        assert seq.seq == ["ACCCACACACCCA", "GGCTTGGAGGAGA"][i]


def test__region_to_seq(genome="tests/data/small_genome.fa.gz", region="chrI:10-20"):
    g = genomepy.Genome(genome)

    # extract sequences marked in track from small_genome.fa.gz
    seq = g._region_to_seq(region=region, extend_up=0, extend_down=0)
    assert seq == "CCCACACACC"

    # extend by varying amounts
    seq = g._region_to_seq(region=region, extend_up=1, extend_down=2)
    # original:    "CCCACACACC"
    assert seq == "ACCCACACACCCA"


def test__regions_to_seqs(
    genome="tests/data/small_genome.fa.gz", track="tests/data/regions.txt"
):
    g = genomepy.Genome(genome)

    # extract sequences marked in regions.bed from small_genome.fa.gz
    seqs = g._regions_to_seqs(track=track, extend_up=0, extend_down=0)
    for i, seq in enumerate(seqs):
        assert seq.name == ["chrI:10-20", "chrII:20-30"][i]
        assert seq.seq == ["CCCACACACC", "TCCTCCAAGC"][i]

    # extend by varying amounts
    seqs = g._regions_to_seqs(track=track, extend_up=1, extend_down=2)
    for i, seq in enumerate(seqs):
        assert seq.name == ["chrI:10-20", "chrII:20-30"][i]
        # original:         "CCCACACACC",    "TCCTCCAAGC"
        assert seq.seq == ["ACCCACACACCCA", "CTCCTCCAAGCCC"][i]


def test_get_track_type():
    tracks = [
        (("chr1:10-20", "chr2:10-20"), "interval"),
        (["chr1:10-20", "chr2:10-20"], "interval"),
        ("tests/data/regions.txt", "interval"),
        ("tests/data/regions.bed", "bed"),
    ]

    for track, track_type in tracks:
        result = genomepy.Genome.get_track_type(track)
        assert result == track_type


def test_track2fasta():
    pass


def test_gap_sizes(genome="tests/data/gap.fa"):
    g = genomepy.Genome(genome)
    g.gap_sizes()

    assert isinstance(g._gap_sizes, dict)
    assert list(g._gap_sizes.keys()) == ["chr1", "chr3"]

    os.unlink(genome + ".fai")
    os.unlink(genome + ".sizes")


def test__weighted_selection():
    pass


def test_get_random_sequences():
    pass


def test_delete_index(genome="tests/data/small_genome.fa.gz"):
    index_file = genome + ".fai"
    os.unlink(index_file)
