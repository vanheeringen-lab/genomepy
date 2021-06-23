import pandas as pd
import pytest

import genomepy.files
import genomepy.utils


def test_create(provider):
    u = provider.create("url")
    # inheritance
    assert u.name == "URL"
    assert hasattr(u, "list_providers")
    # create method overwritten
    assert u.create == NotImplementedError

    with pytest.raises(ValueError):
        provider.create("error")


def test_list_providers(provider):
    assert provider.list_providers() == list(genomepy.provider.PROVIDERS.keys())


def test_online_providers(provider):
    # provider given
    providers = list(provider.online_providers("url"))
    assert len(providers) == 1
    assert providers[0].name == "URL"

    # not given
    providers = list(provider.online_providers())
    assert len(providers) == 4


def test_search_all(provider):
    # unrecognized provider/genome will cause an exception or stopiteration respectively
    # case insensitive description search
    search = list(provider.search_all("xEnOpUs TrOpIcAlIs", "ensembl"))
    metadata = search[0]

    # case insensitive assembly name search
    search = list(provider.search_all("XeNoPuS_tRoPiCaLiS_v9.1", "ensembl"))
    metadata2 = search[0]

    assert metadata == metadata2
    assert isinstance(metadata, list)
    assert "Xenopus_tropicalis_v9.1" in str(metadata[0])
    assert "Ensembl" in str(metadata[1])
    assert "GCA_000004195" in str(metadata[2])
    assert "8364" in str(metadata[3])


def test_download_assembly_report(provider):
    report = provider.download_assembly_report("GCA_000004335")
    assert isinstance(report, pd.DataFrame)
    assert list(report.columns) == genomepy.provider.ASM_FORMAT


def test_map_location(provider):
    readme = "tests/data/ailMel1/README.txt"
    metadata = {
        "name": "ailMel1",
        "provider": "UCSC",
        "assembly_accession": "GCA_000004335.1",
    }
    genomepy.utils.mkdir_p("tests/data/ailMel1")
    genomepy.files.update_readme(readme, metadata)
    genomes_dir = "tests/data"
    mapping = provider.map_locations(
        frm="ailMel1", to="ensembl", genomes_dir=genomes_dir
    )
    assert isinstance(mapping, pd.DataFrame)
    genomepy.utils.rm_rf("tests/data/ailMel1")
