import os

import pytest

import genomepy.files
import genomepy.utils


# to ignore file changes
# git update-index --assume-unchanged tests/data/small_genome.fa.gz
# to recheck file changes
# git update-index --no-assume-unchanged tests/data/small_genome.fa.gz
def test_genome__init__(small_genome):
    # no fasta file
    with pytest.raises(FileNotFoundError):
        genomepy.Genome("empty", "tests/data/genome")

    # genome dir not found
    with pytest.raises(FileNotFoundError):
        genomepy.Genome("unknown", "unknown")

    readme = "tests/data/README.txt"
    if os.path.exists(readme):
        os.unlink(readme)

    assert small_genome.genomes_dir == genomepy.utils.get_genomes_dir(None, False)
    assert small_genome.name == "small_genome"
    assert small_genome.filename == os.path.abspath("tests/data/small_genome.fa.gz")
    assert small_genome.genome_dir == os.path.dirname(small_genome.filename)
    assert os.path.exists(small_genome.index_file)
    assert os.path.exists(small_genome.sizes_file)
    assert os.path.exists(small_genome.gaps_file)
    assert isinstance(small_genome.sizes, dict)
    assert isinstance(small_genome.gaps, dict)
    assert small_genome.annotation_gtf_file is None
    assert small_genome.annotation_bed_file is None
    assert small_genome.tax_id == small_genome.assembly_accession == "na"
    assert isinstance(small_genome.plugin, dict)


def test__parse_name(small_genome):
    # name
    name = small_genome._parse_name("test")
    assert name == "test"

    # file
    name = small_genome._parse_name("/home/genomepy/genomes/test2.fa")
    assert name == "test2"

    # url
    name = small_genome._parse_name(
        "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XT9_1.fa.gz"
    )
    assert name == "XT9_1"


def test__parse_filename(small_genome):
    gpath = "tests/data/small_genome.fa.gz"

    # file path
    filename = small_genome._parse_filename(gpath)
    assert filename == os.path.abspath(gpath)

    # folder path
    filename = small_genome._parse_filename(os.path.dirname(gpath))
    assert filename == os.path.abspath(gpath)

    # name of genome in genomes_dir
    os.mkdir("tests/data/small_genome")
    with open("tests/data/small_genome/small_genome.fa.gz", "w") as fa:
        fa.write("test")
    small_genome.genomes_dir = "tests/data/"
    filename = small_genome._parse_filename(os.path.basename(gpath))
    assert filename == "tests/data/small_genome/small_genome.fa.gz"
    genomepy.utils.rm_rf("tests/data/small_genome")

    # genome not found
    with pytest.raises(FileNotFoundError):
        small_genome._parse_filename("does not exist")


def test__check_annotation_file(small_genome):
    # does not exist
    gtf = small_genome._check_annotation_file("gtf")
    assert gtf is None

    # does exist
    path = "tests/data/small_genome.annotation.test.gz"
    with open(path, "w") as fa:
        fa.write("test")
    test = small_genome._check_annotation_file("test")
    assert test == os.path.abspath(path)
    os.unlink(path)


def test__bed_to_seqs(small_genome, track="tests/data/regions.bed"):
    # extract sequences marked in regions.bed from small_genome.fa.gz
    seqs = small_genome._bed_to_seqs(
        track=track, stranded=False, extend_up=0, extend_down=0
    )
    for i, seq in enumerate(seqs):
        assert seq.name == ["chrI:10-20 gene_a", "chrII:20-30 gene_b"][i]
        assert seq.seq == ["CCCACACACC", "TCCTCCAAGC"][i]

    # second sequence is on the negative strand
    seqs = small_genome._bed_to_seqs(
        track=track, stranded=True, extend_up=0, extend_down=0
    )
    for i, seq in enumerate(seqs):
        assert seq.name == ["chrI:10-20 gene_a", "chrII:20-30 gene_b"][i]
        # original:        "CCCACACACC", "TCCTCCAAGC"
        assert seq.seq == ["CCCACACACC", "GCTTGGAGGA"][i]

    # extend by varying amounts
    seqs = small_genome._bed_to_seqs(
        track=track, stranded=True, extend_up=1, extend_down=2
    )
    for i, seq in enumerate(seqs):
        assert seq.name == ["chrI:10-20 gene_a", "chrII:20-30 gene_b"][i]
        # original:         "CCCACACACC",    "GCTTGGAGGA"
        assert seq.seq == ["ACCCACACACCCA", "GGCTTGGAGGAGA"][i]


def test__region_to_seq(small_genome, region="chrI:10-20"):
    # extract sequences marked in track from small_genome.fa.gz
    seq = small_genome._region_to_seq(region=region, extend_up=0, extend_down=0)
    assert seq == "CCCACACACC"

    # extend by varying amounts
    seq = small_genome._region_to_seq(region=region, extend_up=1, extend_down=2)
    # original:    "CCCACACACC"
    assert seq == "ACCCACACACCCA"


def test__regions_to_seqs(small_genome, track="tests/data/regions.txt"):
    # extract sequences marked in regions.bed from small_genome.fa.gz
    seqs = small_genome._regions_to_seqs(track=track, extend_up=0, extend_down=0)
    for i, seq in enumerate(seqs):
        assert seq.name == ["chrI:10-20", "chrII:20-30"][i]
        assert seq.seq == ["CCCACACACC", "TCCTCCAAGC"][i]

    # extend by varying amounts
    seqs = small_genome._regions_to_seqs(track=track, extend_up=1, extend_down=2)
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
        ("tests/data/regions2.bed", "bed"),
    ]
    for track, track_type in tracks:
        result = genomepy.Genome.get_track_type(track)
        assert result == track_type


def test_track2fasta(small_genome):
    tracks = [
        ("tests/data/regions.txt", "interval"),
        ("tests/data/regions.bed", "bed"),
    ]
    for i, track in enumerate(tracks):
        seq = small_genome.track2fasta(
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


def test_sizes(gap_genome):
    assert list(gap_genome.sizes.keys()) == ["chr1", "chr2", "chr3"]
    assert all(
        isinstance(gap_genome.sizes[chrom], int) for chrom in gap_genome.sizes.keys()
    )
    assert gap_genome.sizes["chr1"] == 28

    # does not overwrite user-set sizes
    gap_genome.sizes = {"asd": 1}
    assert gap_genome.sizes == {"asd": 1}

    # repopulates empty dicts
    gap_genome.sizes = {}
    assert list(gap_genome.sizes.keys()) == ["chr1", "chr2", "chr3"]


def test_gaps(gap_genome):
    assert list(gap_genome.gaps.keys()) == ["chr1", "chr3"]

    # does not overwrite user-set gaps
    gap_genome.gaps = {"asd": 1}
    assert gap_genome.gaps == {"asd": 1}

    # repopulates empty dicts
    gap_genome.gaps = {}
    assert list(gap_genome.gaps.keys()) == ["chr1", "chr3"]


def test__weighted_selection(n=2):
    tuples = [(1, "one"), (2, "two"), (3, "three")]
    ws = genomepy.Genome._weighted_selection(tuples, n)

    assert len(ws) == n
    assert isinstance(ws[0], str)
    assert tuples[0][1] in ws or tuples[1][1] in ws or tuples[2][1] in ws


def test_get_random_sequences(small_genome):
    n = 2
    length = 200  # default
    chroms = ["chrI", "chrII"]
    max_n = 0.1  # default
    rs = small_genome.get_random_sequences(
        n=n, length=length, chroms=chroms, max_n=max_n
    )

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
    rs = small_genome.get_random_sequences(n=1, chroms=chroms, outtype="string")
    assert (
        str(small_genome.track2fasta(rs[0])[0].seq).upper().count("N") <= length * max_n
    )