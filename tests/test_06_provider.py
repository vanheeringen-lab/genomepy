import os
from tempfile import TemporaryDirectory

import pandas as pd
import pytest

import genomepy
import genomepy.utils
from tests import linux, travis
from tests.conftest import validate_annot


def test_provider_status(provider):
    provider.provider_status("https://www.google.com")

    with pytest.raises(ConnectionError):
        provider.provider_status("https://www.thiswebsiteisoffline.nl/")


def test_providerbase__init__(provider):
    assert list(provider._providers) == ["ensembl", "ucsc", "ncbi", "url"]
    assert provider.name is None
    assert isinstance(provider.genomes, dict)
    assert isinstance(provider.accession_fields, list)


def test_create(provider):
    provider.create("url")

    with pytest.raises(ValueError):
        provider.create("error")


def test_register_provider_and_list_providers(provider):
    assert isinstance(provider._providers, dict)

    providers = ["ensembl", "ucsc", "ncbi", "url"]
    for p in providers:
        assert p in list(provider.list_providers())


def test__genome_info_tuple(provider):
    with pytest.raises(NotImplementedError):
        provider._genome_info_tuple(None)


def test_list_available_genomes(provider):
    provider.list_available_genomes()


def test_online_providers(provider):
    # provider given
    providers = list(provider.online_providers("url"))
    assert len(providers) == 1
    assert providers[0].name == "URL"

    # not given
    providers = list(provider.online_providers())
    assert (
        len(providers) == 4
    ), "Looks like a provider is offline. Not your fault (this time)!"


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
    assert "8364" in str(metadata[4])


def test_check_name(ucsc):
    ucsc.check_name("ailMel1")
    with pytest.raises(genomepy.exceptions.GenomeDownloadError):
        ucsc.check_name("not_a_real_genome")


def get_genome_download_link(provider):
    with pytest.raises(NotImplementedError):
        provider.get_genome_download_link()


def test_genome_taxid(ucsc):
    taxid = ucsc.genome_taxid(ucsc.genomes["ailMel1"])
    assert taxid == 9646


def test_assembly_accession(ucsc):
    accession = ucsc.assembly_accession(ucsc.genomes["ailMel1"])
    assert "000004335" in accession


def test_download_assembly_report(provider):
    report = provider.download_assembly_report("GCA_000004335")
    assert isinstance(report, pd.DataFrame)
    assert list(report.columns) == genomepy.provider.ASM_FORMAT


@pytest.mark.skipif(not travis, reason="slow")
def test_download_genome(
    ucsc,
    name="sacCer3",
    localname="my_genome",
    mask="soft",
):
    out_dir = os.getcwd()

    with TemporaryDirectory(dir=out_dir) as tmpdir:
        ucsc.download_genome(
            name,
            genomes_dir=tmpdir,
            localname=localname,
            mask=mask,
        )

        genome = os.path.join(tmpdir, localname, localname + ".fa")
        assert os.path.exists(genome)

        readme = os.path.join(os.path.dirname(genome), "README.txt")
        with open(readme) as f:
            metadata = {}
            for line in f.readlines():
                vals = line.strip().split(":")
                metadata[vals[0].strip()] = (":".join(vals[1:])).strip()

        assert metadata["name"] == localname
        assert metadata["mask"] == mask


def test_map_location(provider):
    readme = "tests/data/ailMel1/README.txt"
    metadata = {
        "name": "ailMel1",
        "provider": "UCSC",
        "assembly_accession": "GCA_000004335.1",
    }
    genomepy.utils.mkdir_p("tests/data/ailMel1")
    genomepy.utils.update_readme(readme, metadata)
    genomes_dir = "tests/data"
    mapping = provider.map_locations(
        frm="ailMel1", to="ensembl", genomes_dir=genomes_dir
    )
    assert isinstance(mapping, pd.DataFrame)
    genomepy.utils.rm_rf("tests/data/ailMel1")


def test_get_annotation_download_link(provider):
    with pytest.raises(NotImplementedError):
        provider.get_annotation_download_link(None)


@pytest.mark.skipif(not travis, reason="slow")
@pytest.mark.skipif(travis and linux, reason="FTP does not work on Travis-Linux")
def test_download_and_generate_annotation(provider):
    out_dir = os.getcwd()
    localname = "my_annot"

    annot_url = "https://www.google.com"
    with pytest.raises(TypeError), TemporaryDirectory(dir=out_dir) as tmpdir:
        provider.download_and_generate_annotation(
            genomes_dir=tmpdir, annot_url=annot_url, localname=localname
        )

    # FTP used here to test FTP file downloading
    annot_url = (
        "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/"
        + "GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.gff.gz"
    )
    with TemporaryDirectory(dir=out_dir) as tmpdir:
        provider.download_and_generate_annotation(
            genomes_dir=tmpdir, annot_url=annot_url, localname=localname
        )

        fname = os.path.join(tmpdir, localname, localname + ".annotation.gtf")
        validate_annot(fname, "gtf")

        fname = os.path.join(tmpdir, localname, localname + ".annotation.bed")
        validate_annot(fname, "bed")


def test_attempt_and_report(provider, caplog):
    out_dir = os.getcwd()
    localname = "my_annot"
    name = "test"

    provider.attempt_and_report(name, localname, None, None)
    assert f"Could not download gene annotation for {name} from" in caplog.text

    annot_url = "https://www.google.com"
    with pytest.raises(genomepy.exceptions.GenomeDownloadError), TemporaryDirectory(
        dir=out_dir
    ) as tmpdir:
        provider.attempt_and_report(name, localname, annot_url, tmpdir)

    assert f"Downloading annotation from None. Target URL: {annot_url}" in caplog.text


@pytest.mark.skipif(not travis, reason="slow")
def test_download_annotation(ucsc):
    out_dir = os.getcwd()
    localname = "my_annot"

    name = "xenTro2"

    annot_urls = [
        "http://hgdownload.cse.ucsc.edu/goldenPath/xenTro2/database/ensGene.txt.gz",
        "http://hgdownload.cse.ucsc.edu/goldenPath/xenTro2/database/refGene.txt.gz",
    ]
    with TemporaryDirectory(dir=out_dir) as tmpdir:
        ucsc.download_annotation(name=name, genomes_dir=tmpdir, localname=localname)

        # check download_and_generate_annotation output
        fname = os.path.join(tmpdir, localname, localname + ".annotation.gtf")
        validate_annot(fname, "gtf")

        fname = os.path.join(tmpdir, localname, localname + ".annotation.bed")
        validate_annot(fname, "bed")

        # check attempt_download_and_report_back output
        readme = os.path.join(tmpdir, localname, "README.txt")
        metadata, lines = genomepy.utils.read_readme(readme)
        assert metadata["annotation url"] in annot_urls


def test__search_taxids(ucsc):
    assert not ucsc._search_taxids(ucsc.genomes["ailMel1"], "not_an_id")
    assert ucsc._search_taxids(ucsc.genomes["ailMel1"], "9646")


def test__search_descriptions(ucsc):
    assert "scientificName" in ucsc.description_fields
    assert ucsc.genomes["ailMel1"]["scientificName"] == "Ailuropoda melanoleuca"
    desc = genomepy.utils.safe("Ailuropoda melanoleuca").lower()
    assert ucsc._search_descriptions(ucsc.genomes["ailMel1"], desc)
    assert not ucsc._search_descriptions(ucsc.genomes["ailMel1"], "not_in_description")


def test_search(ucsc):
    for method in ["ailMel1", "9646", "Ailuropoda melanoleuca"]:
        genomes = ucsc.search(method)
        for genome in genomes:
            assert genome[0] == "ailMel1"
