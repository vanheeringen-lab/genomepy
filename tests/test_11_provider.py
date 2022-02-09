import pandas as pd
import pytest

import genomepy.files
import genomepy.utils
from tests import linux, travis


def test_create():
    u = genomepy.providers.create("url")
    assert u.name == "URL"
    assert hasattr(u, "download_assembly_report")

    with pytest.raises(ValueError):
        genomepy.providers.create("error")


def test_list_providers():
    assert len(genomepy.providers.list_providers()) == len(genomepy.providers.PROVIDERS)


def test_list_online_providers():
    _all = genomepy.providers.list_providers()
    online = genomepy.providers.list_online_providers()
    if linux and travis:
        assert online == _all[1:]  # GENCODE is on FTP
    else:
        assert online == _all


def test_online_providers():
    # provider given
    providers = list(genomepy.providers.online_providers("url"))
    assert len(providers) == 1
    assert providers[0].name == "URL"

    # not given
    providers = list(genomepy.providers.online_providers())
    if linux and travis:
        assert len(providers) == 5  # GENCODE is on FTP
    else:
        assert len(providers) == 6


def test_search():
    # unrecognized provider/genome will cause an exception or stopiteration respectively
    # case insensitive description search
    search = list(genomepy.providers.search("xEnOpUs TrOpIcAlIs", "ensembl"))
    metadata = search[0]

    # case insensitive assembly name search
    search = list(genomepy.providers.search("XeNoPuS_tRoPiCaLiS_v9.1", "ensembl"))
    metadata2 = search[0]

    assert metadata == metadata2
    assert isinstance(metadata, list)
    assert "Xenopus_tropicalis_v9.1" in str(metadata[0])
    assert "Ensembl" in str(metadata[1])
    assert "GCA_000004195" in str(metadata[2])
    assert "8364" in str(metadata[3])


def test_map_location():
    mapping = genomepy.providers.map_locations(
        frm="sacCer3", to="ensembl", genomes_dir="tests/data"
    )
    assert isinstance(mapping, pd.DataFrame)
    assert mapping.index.name == "ucsc_name"
    assert mapping.columns.to_list() == ["ensembl_name"]


def test_provider(provider):
    assert hasattr(provider, "list")
    assert hasattr(provider, "online")
    assert hasattr(provider, "create")
