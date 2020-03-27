import genomepy
import os


# to ignore file changes
# git update-index --assume-unchanged tests/data/small_genome.fa.gz
# to recheck file changes
# git update-index --no-assume-unchanged tests/data/small_genome.fa.gz
def test_genome_index_generation(genome="tests/data/small_genome.fa.gz"):
    index_file = genome + ".fai"
    if os.path.exists(index_file):
        os.unlink(index_file)

    # initialize the class (creates the index file
    genomepy.Genome(genome)
    assert os.path.exists(index_file)
    os.unlink(index_file)


def test__bed_to_seqs():
    pass


def test__region_to_seqs():
    pass


def test_get_track_type():
    pass


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
