from tempfile import NamedTemporaryFile
from genomepy.utils import generate_gap_bed, generate_fa_sizes


def test_gaps():
    infa = "tests/data/gap.fa"
    outbed = "tests/data/gap.bed"

    tmp = NamedTemporaryFile().name
    generate_gap_bed(infa, tmp)

    result = open(tmp).read()
    expect = open(outbed).read()

    assert result == expect


def test_sizes():
    infa = "tests/data/gap.fa"

    tmp = NamedTemporaryFile().name
    generate_fa_sizes(infa, tmp)

    result = open(tmp).read()

    assert result == "chr1\t28\nchr2\t45\nchr3\t15\n"
