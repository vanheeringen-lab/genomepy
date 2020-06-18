import genomepy
import os
import pytest
import shutil

from stat import S_IREAD, S_IRGRP, S_IROTH
from genomepy.provider import ProviderBase


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

    readme = "tests/data/README.txt"
    if os.path.exists(readme):
        os.unlink(readme)

    g = genomepy.Genome(genome)
    assert g.genomes_dir == genomepy.utils.get_genomes_dir(None, False)
    assert g.name == "small_genome"
    assert g.filename == os.path.abspath(genome)
    assert g.genome_dir == os.path.dirname(g.filename)
    assert os.path.exists(g.index_file)
    assert os.path.exists(g.sizes_file)
    assert os.path.exists(g.gaps_file)
    assert isinstance(g.sizes, dict)
    assert isinstance(g.gaps, dict)
    assert g.annotation_gtf_file is None
    assert g.annotation_bed_file is None
    assert g.tax_id == g.assembly_accession == "na"
    assert isinstance(g.plugin, dict)


def test__parse_name(genome="tests/data/small_genome.fa.gz"):
    g = genomepy.Genome(genome)  # unimportant

    # name
    name = g._parse_name("test")
    assert name == "test"

    # file
    name = g._parse_name("/home/genomepy/genomes/test2.fa")
    assert name == "test2"

    # url
    name = g._parse_name("http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XT9_1.fa.gz")
    assert name == "XT9_1"


def test__parse_filename(genome="tests/data/small_genome.fa.gz"):
    g = genomepy.Genome(genome)  # unimportant

    # file path
    filename = g._parse_filename(genome)
    assert filename == os.path.abspath(genome)

    # folder path
    filename = g._parse_filename(os.path.dirname(genome))
    assert filename == os.path.abspath(genome)

    # name of genome in genomes_dir
    os.mkdir("tests/data/small_genome")
    with open("tests/data/small_genome/small_genome.fa.gz", "w") as fa:
        fa.write("test")
    g.genomes_dir = "tests/data/"
    filename = g._parse_filename(os.path.basename(genome))
    assert filename == "tests/data/small_genome/small_genome.fa.gz"
    shutil.rmtree("tests/data/small_genome")

    # genome not found
    with pytest.raises(FileNotFoundError):
        g._parse_filename("does not exist")


def test_check_annotation_file(genome="tests/data/small_genome.fa.gz"):
    g = genomepy.Genome(genome)

    # does not exist
    gtf = g.check_annotation_file("gtf")
    assert gtf is None

    # does exist
    path = "tests/data/small_genome.annotation.test.gz"
    with open(path, "w") as fa:
        fa.write("test")
    test = g.check_annotation_file("test")
    assert test == os.path.abspath(path)
    os.unlink(path)


def test__update_provider(genome="tests/data/small_genome.fa.gz"):
    g = genomepy.Genome(genome)

    # can't parse url
    metadata = {}
    g._update_provider(metadata)
    assert metadata.get("provider") == "Unknown"

    # can parse url
    metadata = {
        "genome url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/465/"
        "GCF_000146465.1_ASM14646v1/GCF_000146465.1_ASM14646v1_genomic.fna.gz"
    }
    g._update_provider(metadata)
    assert metadata.get("provider") == "NCBI"


def test__update_tax_id(genome="tests/data/small_genome.fa.gz"):
    g = genomepy.Genome(genome)

    # genome not found
    metadata = {}
    g._update_tax_id(metadata)
    assert metadata["tax_id"] == "na"

    # genome found
    metadata = {}
    provider = ProviderBase.create("NCBI")
    genome = provider.genomes.get("ASM14646v1")

    g._update_tax_id(metadata, provider, genome)
    assert metadata["tax_id"] == "58839"


def test__update_assembly_accession(genome="tests/data/small_genome.fa.gz"):
    g = genomepy.Genome(genome)

    # genome not found
    metadata = {}
    g._update_assembly_accession(metadata)
    assert metadata["assembly_accession"] == "na"

    # genome found
    metadata = {}
    provider = ProviderBase.create("NCBI")
    genome = provider.genomes.get("ASM14646v1")

    g._update_assembly_accession(metadata, provider, genome)
    assert metadata["assembly_accession"] == "GCA_000146465.1"


def test__update_metadata(genome="tests/data/small_genome.fa.gz"):
    g = genomepy.Genome(genome)

    metadata = {"provider": "NCBI", "original name": "ASM14646v1"}
    g._update_metadata(metadata)
    assert metadata["tax_id"] == "58839"
    assert metadata["assembly_accession"] == "GCA_000146465.1"


def test__read_metadata(genome="tests/data/small_genome.fa.gz"):
    g = genomepy.Genome(genome)

    # no readme found
    readme = g.readme_file
    if os.path.exists(readme):
        os.unlink(readme)
    metadata = g._read_metadata()
    assert metadata["provider"] == "na"

    # no overwrites to metadata
    with open(readme, "w") as f:
        f.writelines("provider: not_really_NCBI\n")
        f.writelines("tax_id: not_really_58839\n")
        f.writelines("assembly_accession: not_really_GCA_000146465.1\n")
    metadata = g._read_metadata()
    assert metadata["provider"] == "not_really_NCBI"

    # updates to metadata dict and file
    with open(readme, "w") as f:
        f.writelines(
            "genome url: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/"
            "146/465/GCF_000146465.1_ASM14646v1/"
            "GCF_000146465.1_ASM14646v1_genomic.fna.gz\n"
        )
        f.writelines("tax_id: not_really_58839\n")
        f.writelines("assembly_accession: not_really_GCA_000146465.1\n")
    metadata1 = g._read_metadata()
    assert metadata1["provider"] == "NCBI"
    metadata2, _ = genomepy.utils.read_readme(readme)
    assert metadata2["provider"] == "NCBI"

    # no writing permission to file
    with open(readme, "w") as f:
        f.writelines(
            "genome url: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/"
            "146/465/GCF_000146465.1_ASM14646v1/"
            "GCF_000146465.1_ASM14646v1_genomic.fna.gz\n"
        )
    os.chmod(readme, S_IREAD | S_IRGRP | S_IROTH)
    metadata1 = g._read_metadata()
    assert metadata1["provider"] == "na"
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
    assert list(g.sizes.keys()) == ["chr1", "chr2", "chr3"]
    assert all(isinstance(g.sizes[chrom], int) for chrom in g.sizes.keys())
    assert g.sizes["chr1"] == 28

    # does not overwrite user-set sizes
    g.sizes = {"asd": 1}
    assert g.sizes == {"asd": 1}

    # repopulates empty dicts
    g.sizes = {}
    assert list(g.sizes.keys()) == ["chr1", "chr2", "chr3"]


def test_gaps(genome="tests/data/gap.fa"):
    g = genomepy.Genome(genome)
    assert list(g.gaps.keys()) == ["chr1", "chr3"]

    # does not overwrite user-set gaps
    g.gaps = {"asd": 1}
    assert g.gaps == {"asd": 1}

    # repopulates empty dicts
    g.gaps = {}
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
