import genomepy
import os
import pytest

from platform import system


linux = system() == "Linux"
travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"


@pytest.fixture(scope="session")
def ce11_genome(tmpdir_factory):
    tmpdirname = tmpdir_factory.mktemp("genomes")
    genomepy.functions.install_genome(
        name="ce11", provider="UCSC", genomes_dir=tmpdirname, annotation=True,
    )

    g = genomepy.Genome("ce11", genomes_dir=tmpdirname)
    return g


@pytest.mark.skipif(travis or not linux, reason="slow")
def test_gene_annotation(ce11_genome):
    df = ce11_genome.gene_annotation()
    assert "chrV" in df["chrom"].values


@pytest.mark.skipif(travis or not linux, reason="slow")
def test_ensembl_info(ce11_genome):
    info = ce11_genome.ensembl_genome_info()
    assert info == ("WBcel235", "GCA_000002985.3", "6239")
