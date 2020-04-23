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

    g = genomepy.Genome(genome)
    assert os.path.exists(g.index_file)
    assert not g.annotation_gtf_file
    assert not g.annotation_bed_file
    assert os.path.exists(g.sizes_file)
    assert os.path.exists(g.gaps_file)
    assert g.tax_id == g.assembly_accession == "na"


def test__read_metadata(capsys, genome="tests/data/small_genome.fa.gz"):
    # create blank README.txt
    g = genomepy.Genome(genome)
    readme = g.readme_file
    if os.path.exists(readme):
        os.unlink(readme)
    metadata, lines = genomepy.utils.read_readme(readme)

    assert metadata["provider"] == "na"

    # no changes to metadata
    with open(readme, "w") as f:
        f.writelines("provider: NCBI\n")
        f.writelines("original name: ASM14646v1\n")
        f.writelines("tax_id: 58839\n")
        f.writelines("assembly_accession: GCA_000146465.1\n")
    genomepy.Genome(genome)

    metadata, lines = genomepy.utils.read_readme(readme)
    assert metadata["provider"] == "NCBI"
    assert metadata["original name"] == "ASM14646v1"
    assert metadata["tax_id"] == "58839"
    assert metadata["assembly_accession"] == "GCA_000146465.1"
    os.unlink(readme)


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


def test_track2fasta(genome="tests/data/small_genome.fa.gz"):
    tracks = [
        ("tests/data/regions.txt", "interval"),
        ("tests/data/regions.bed", "bed"),
    ]
    g = genomepy.Genome(genome)

    for i, track in enumerate(tracks):
        seq = g.track2fasta(
            track=track[0],
            fastafile=None,
            stranded=False,
            extend_up=i,
            extend_down=i + 1,
        )

        # default sequence:       CCCACACACC
        if i == 0:  # extend up +0, down -1
            assert seq[0].seq == "CCCACACACCC"
            assert seq[1].seq == "TCCTCCAAGCC"
        else:  # extend up +1, down -4
            assert seq[0].seq == "ACCCACACACCCA"
            assert seq[1].seq == "CTCCTCCAAGCCC"


def test_sizes(genome="tests/data/gap.fa"):
    g = genomepy.Genome(genome)
    g.contig_sizes()
    assert list(g.sizes.keys()) == ["chr1", "chr2", "chr3"]


def test_gaps(genome="tests/data/gap.fa"):
    g = genomepy.Genome(genome)
    g.gap_sizes()
    assert list(g.gaps.keys()) == ["chr1", "chr3"]


def test__weighted_selection(n=2):
    tuples = [(1, "one"), (2, "two"), (3, "three")]
    ws = genomepy.Genome._weighted_selection(tuples, n)

    assert len(ws) == n
    assert isinstance(ws[0], str)
    assert tuples[0][1] in ws or tuples[1][1] in ws or tuples[2][1] in ws


def test_get_random_sequences(genome="tests/data/small_genome.fa.gz"):
    g = genomepy.Genome(genome)
    n = 2
    length = 200  # default
    chroms = ["chrI", "chrII"]
    max_n = 0.1  # default
    rs = g.get_random_sequences(n=n, length=length, chroms=chroms, max_n=max_n)

    # check that the output has the right length, content, types, and sequence length
    assert len(rs) == n
    for i in range(n):
        assert rs[i][0] in chroms
        assert (
            isinstance(rs[i][0], str)
            and isinstance(rs[i][1], int)
            and isinstance(rs[i][2], int)
        )
        assert rs[i][2] - rs[i][1] == length

    # check that the max Ns are lower than the expected cutoff
    rs = g.get_random_sequences(n=1, chroms=chroms, outtype="string")
    assert str(g.track2fasta(rs[0])[0].seq).upper().count("N") <= length * max_n


def test_delete_test_files():
    for genome in [
        "tests/data/small_genome.",
        "tests/data/gap.",
    ]:
        for ext in ["fa.fai", "fa.sizes", "gaps.bed", "fa.gz.fai", "fa.gz.sizes"]:
            file = genome + ext
            if os.path.exists(file):
                os.unlink(file)
