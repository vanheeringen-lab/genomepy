import genomepy
import pytest


@pytest.fixture(
    scope="module",
    params=[
        ["GCA_000004195.3", "NCBI", "Xenopus_tropicalis_v9.1"],
        ["GCA_000004195.3", "Ensembl", "Xenopus_tropicalis_v9.1"],
        ["GCA_000004195.3", "UCSC", "xenTro9"],
    ],
)
def accession_data(request):
    return request.param


@pytest.fixture(
    scope="module",
    params=[
        ["8364", "NCBI", "Xenopus_tropicalis_v9.1"],
        ["8364", "Ensembl", "Xenopus_tropicalis_v9.1"],
        ["8364", "UCSC", "xenTro9"],
    ],
)
def taxid_data(request):
    return request.param


def test_assembly_accession(accession_data):
    """Test retrieving assembly accession for a given genome genome name"""
    acc, provider, name = accession_data

    p = genomepy.provider.ProviderBase.create(provider)
    assert acc == p.assembly_accession(name)


def test_taxid(taxid_data):
    """Test retrieving taxonomy_id for a given genome genome name"""
    taxid, provider, name = taxid_data

    p = genomepy.provider.ProviderBase.create(provider)
    assert taxid == str(p.genome_taxid(name))
