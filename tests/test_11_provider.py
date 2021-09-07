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
        assert len(providers) == 4  # GENCODE is on FTP
    else:
        assert len(providers) == 5


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


def test__best_search_result():
    # ncbi_human = list(genomepy.search("GCA_000001405", provider="NCBI"))
    ncbi_human = [
        [
            "GRCh38",
            "NCBI",
            "GCF_000001405.26",
            9606,
            True,
            "Homo sapiens",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p1",
            "NCBI",
            "GCF_000001405.27",
            9606,
            True,
            "Homo sapiens",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p2",
            "NCBI",
            "GCF_000001405.28",
            9606,
            True,
            "Homo sapiens",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p3",
            "NCBI",
            "GCF_000001405.29",
            9606,
            True,
            "Homo sapiens",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p4",
            "NCBI",
            "GCF_000001405.30",
            9606,
            True,
            "Homo sapiens",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p5",
            "NCBI",
            "GCF_000001405.31",
            9606,
            True,
            "Homo sapiens",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p6",
            "NCBI",
            "GCF_000001405.32",
            9606,
            True,
            "Homo sapiens",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p7",
            "NCBI",
            "GCF_000001405.33",
            9606,
            True,
            "Homo sapiens",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p8",
            "NCBI",
            "GCF_000001405.34",
            9606,
            True,
            "Homo sapiens",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p9",
            "NCBI",
            "GCF_000001405.35",
            9606,
            True,
            "Homo sapiens",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p10",
            "NCBI",
            "GCF_000001405.36",
            9606,
            True,
            "Homo sapiens",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p11",
            "NCBI",
            "GCF_000001405.37",
            9606,
            True,
            "Homo sapiens",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p12",
            "NCBI",
            "GCF_000001405.38",
            9606,
            True,
            "Homo sapiens",
            "Genome Reference Consortium",
        ],
        [
            "GRCh37.p13",
            "NCBI",
            "GCF_000001405.25",
            9606,
            True,
            "Homo sapiens",
            "Genome Reference Consortium",
        ],
        [
            "GRCh38.p13",
            "NCBI",
            "GCF_000001405.39",
            9606,
            True,
            "Homo sapiens",
            "Genome Reference Consortium",
        ],
    ]

    # list of search results
    best_result = genomepy.providers._best_search_result("GCA_000001405.39", ncbi_human)
    assert best_result[2] == "GCF_000001405.39"

    # zero search results
    assert genomepy.providers._best_search_result("GCA_000001405.39", []) is None

    # one search results
    best_result = genomepy.providers._best_search_result(
        "GCA_000001405.26", [ncbi_human[0]]
    )
    assert best_result[2] == ncbi_human[0][2]


def test_download_assembly_report():
    assembly_report = "tests/data/sacCer3/assembly_report.txt"
    genomepy.providers.download_assembly_report("GCA_000146045", assembly_report)
    report = pd.read_csv(assembly_report, sep="\t", comment="#")

    assert isinstance(report, pd.DataFrame)
    assert list(report.columns) == genomepy.providers.ASM_FORMAT


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
