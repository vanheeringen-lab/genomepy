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
    assembly_report = "tests/data/sacCer3/assembly_report.txt"
    provider.download_assembly_report("GCA_000146045", assembly_report)
    report = pd.read_csv(assembly_report, sep="\t", comment="#")

    assert isinstance(report, pd.DataFrame)
    assert list(report.columns) == genomepy.provider.ASM_FORMAT


def test_map_location(provider):
    mapping = provider.map_locations(
        frm="sacCer3", to="ensembl", genomes_dir="tests/data"
    )
    assert isinstance(mapping, pd.DataFrame)
    assert mapping.index.name == "ucsc_name"
    assert mapping.columns.to_list() == ["ensembl_name"]
