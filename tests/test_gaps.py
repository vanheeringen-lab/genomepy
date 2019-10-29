from tempfile import NamedTemporaryFile
from genomepy.utils import generate_gap_bed


def test_gaps():
    infa = "tests/data/gap.fa"
    outbed = "tests/data/gap.bed"

    tmp = NamedTemporaryFile().name
    generate_gap_bed(infa, tmp)

    result = open(tmp).read()
    expect = open(outbed).read()

    assert result == expect
