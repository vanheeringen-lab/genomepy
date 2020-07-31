import genomepy
import gzip
import os
import pytest

from tempfile import TemporaryDirectory
from platform import system

linux = system() == "Linux"
travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"


def validate_gzipped_gtf(fname):
    assert os.path.exists(fname)
    with gzip.open(fname, "r") as f:
        for line in f:
            line = line.decode()
            if line.startswith("#"):
                continue
            vals = line.split("\t")
            assert 9 == len(vals)
            int(vals[3]), int(vals[4])
            break


def validate_gzipped_bed(fname):
    assert os.path.exists(fname)
    with gzip.open(fname, "r") as f:
        for line in f:
            line = line.decode()
            if line.startswith("#"):
                continue
            vals = line.split("\t")
            assert 12 == len(vals)
            int(vals[1]), int(vals[2])
            break


@pytest.fixture(scope="module")
def p():
    return genomepy.provider.ProviderBase()


def test_providerbase__init__(p):
    assert list(p._providers) == ["ensembl", "ucsc", "ncbi", "url"]
    assert p.name is None
    assert isinstance(p.genomes, dict)
    assert isinstance(p.accession_fields, list)


def test_create(p):
    p.create("Ensembl")

    with pytest.raises(ValueError):
        p.create("error")


@pytest.mark.disable_socket
def test_create_offline(p):
    providers = ["ensembl", "ucsc", "ncbi"]
    for provider in providers:
        with pytest.raises(ConnectionError):
            p.create(provider)


def test_register_provider_and_list_providers(p):
    assert isinstance(p._providers, dict)

    providers = ["ensembl", "ucsc", "ncbi", "url"]
    for provider in providers:
        assert provider in list(p.list_providers())


def test__genome_info_tuple(p):
    with pytest.raises(NotImplementedError):
        p._genome_info_tuple(None)


def test_list_available_genomes(p):
    p.list_available_genomes()


def test_check_name(p):
    p = p.create("Ensembl")
    p.check_name("KH")
    with pytest.raises(genomepy.exceptions.GenomeDownloadError):
        p.check_name("not_a_real_genome")


def get_genome_download_link(p):
    with pytest.raises(NotImplementedError):
        p.get_genome_download_link()


def test_genome_taxid(p):
    p = p.create("Ensembl")
    taxid = p.genome_taxid(p.genomes["KH"])
    assert taxid == 7719


def test_assembly_accession(p):
    p = p.create("Ensembl")
    accession = p.assembly_accession(p.genomes["KH"])
    assert accession.startswith("GCA_000224145")


@pytest.mark.skipif(not travis or not linux, reason="slow")
def test_download_genome(
    p,
    name="sacCer3",
    localname="my_genome",
    mask="soft",
    regex="MT",
    invert_match=True,
    bgzip=True,
):
    p = p.create("UCSC")
    out_dir = os.getcwd()

    with TemporaryDirectory(dir=out_dir) as tmpdir:
        p.download_genome(
            name,
            genomes_dir=tmpdir,
            localname=localname,
            mask=mask,
            regex=regex,
            invert_match=invert_match,
            bgzip=bgzip,
        )

        genome = os.path.join(tmpdir, localname, localname + ".fa.gz")
        assert os.path.exists(genome)

        readme = os.path.join(os.path.dirname(genome), "README.txt")
        with open(readme) as f:
            metadata = {}
            for line in f.readlines():
                vals = line.strip().split(":")
                metadata[vals[0].strip()] = (":".join(vals[1:])).strip()

        assert metadata["name"] == localname
        assert metadata["mask"] == mask


def test_get_annotation_download_link(p):
    with pytest.raises(NotImplementedError):
        p.get_annotation_download_link(None)


@pytest.mark.skipif(not travis or not linux, reason="slow")
def test_download_and_generate_annotation(p):
    out_dir = os.getcwd()
    localname = "my_annot"

    annot_url = "https://www.google.com"
    with pytest.raises(TypeError), TemporaryDirectory(dir=out_dir) as tmpdir:
        p.download_and_generate_annotation(
            genomes_dir=tmpdir, annot_url=annot_url, localname=localname
        )

    annot_url = (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/"
        + "GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.gff.gz"
    )
    with TemporaryDirectory(dir=out_dir) as tmpdir:
        p.download_and_generate_annotation(
            genomes_dir=tmpdir, annot_url=annot_url, localname=localname
        )

        fname = os.path.join(tmpdir, localname, localname + ".annotation.gtf.gz")
        validate_gzipped_gtf(fname)

        fname = os.path.join(tmpdir, localname, localname + ".annotation.bed.gz")
        validate_gzipped_bed(fname)


def test_attempt_and_report(p, capsys):
    out_dir = os.getcwd()
    localname = "my_annot"
    name = "test"

    p.attempt_and_report(name, localname, None, None)
    captured = capsys.readouterr().err.strip()
    assert captured.startswith(f"Could not download genome annotation for {name} from")

    annot_url = "https://www.google.com"
    with pytest.raises(genomepy.exceptions.GenomeDownloadError), TemporaryDirectory(
        dir=out_dir
    ) as tmpdir:
        p.attempt_and_report(name, localname, annot_url, tmpdir)

    captured = capsys.readouterr().err.strip()
    assert captured.startswith(f"Downloading annotation from {annot_url}")


@pytest.mark.skipif(not travis or not linux, reason="slow")
def test_download_annotation(p):
    out_dir = os.getcwd()
    localname = "my_annot"

    p = p.create("UCSC")
    name = "xenTro2"

    annot_url = (
        "http://hgdownload.cse.ucsc.edu/goldenPath/xenTro2/database/ensGene.txt.gz"
    )
    with TemporaryDirectory(dir=out_dir) as tmpdir:
        p.download_annotation(name=name, genomes_dir=tmpdir, localname=localname)

        # check download_and_generate_annotation output
        fname = os.path.join(tmpdir, localname, localname + ".annotation.gtf.gz")
        validate_gzipped_gtf(fname)

        fname = os.path.join(tmpdir, localname, localname + ".annotation.bed.gz")
        validate_gzipped_bed(fname)

        # check attempt_download_and_report_back output
        readme = os.path.join(tmpdir, localname, "README.txt")
        metadata, lines = genomepy.utils.read_readme(readme)
        assert metadata["annotation url"] == annot_url


def test__search_taxids(p):
    p = p.create("Ensembl")
    assert not p._search_taxids(p.genomes["KH"], "not_an_id")
    assert p._search_taxids(p.genomes["KH"], "7719")


def test__search_descriptions(p):
    p = p.create("Ensembl")
    assert "scientific_name" in p.description_fields
    assert p.genomes["KH"]["scientific_name"] == "Ciona intestinalis"
    desc = genomepy.utils.safe("Ciona intestinalis").lower()
    assert p._search_descriptions(p.genomes["KH"], desc)
    assert not p._search_descriptions(p.genomes["KH"], "not_in_description")


def test_search(p):
    p = p.create("Ensembl")
    for method in ["KH", "7719", "Ciona intestinalis"]:
        genomes = p.search(method)
        for genome in genomes:
            assert genome[0] == "KH"
