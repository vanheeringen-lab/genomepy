import os
from tempfile import TemporaryDirectory

import pytest

import genomepy
import genomepy.utils
from tests import travis
from tests.conftest import validate_annot


def test_ncbiprovider__init__(url):
    assert url.name == "URL"
    assert url.genomes == {}


def test_genome_taxid(url):
    assert url.genome_taxid({}) == "na"


def test_assembly_accession(url):
    assert url.assembly_accession({}) == "na"


def test_search(url, ucsc):
    result = url.search("???")
    expected = ucsc.search("???")
    assert isinstance(result, type(expected))
    with pytest.raises(StopIteration):
        assert next(result)


def test_get_genome_download_link(url):
    link = url.get_genome_download_link("url")
    assert link == "url"


def test_get_annotation_download_link(url):
    target = "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_Xenbase.gtf"
    link = url.get_annotation_download_link(None, **{"to_annotation": target})
    assert link == target

    with pytest.raises(TypeError):
        bad_url = "bad_url"
        url.get_annotation_download_link(None, **{"to_annotation": bad_url})

    with pytest.raises(TypeError):
        bad_url = "http://good_url.bad_ext"
        url.get_annotation_download_link(None, **{"to_annotation": bad_url})


def test_search_url_for_annotations(url):
    target = "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_genome.fa.gz"
    links = url.search_url_for_annotations(target, "XENTR_9.1")
    expected = [
        "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_Xenbase.gtf",
        "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_GCA.gff3",
        "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_GCF.gff3",
        "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_Xenbase.gff3",
        "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_Xenbase_longest.gff3",
    ]
    assert links == expected

    target = (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/"
        "GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz"
    )
    links = url.search_url_for_annotations(target, "GCF_000027325.1_ASM2732v1")
    expected = [
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/"
        + "GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.gtf.gz",
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/"
        + "GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.gff.gz",
    ]
    assert links == expected

    # no annot file
    with pytest.raises(FileNotFoundError):
        target = (
            "http://ftp.ensembl.org/pub/release-100/fasta/marmota_marmota_marmota/"
            "dna/Marmota_marmota_marmota.marMar2.1.dna.toplevel.fa.gz"
        )
        url.search_url_for_annotations(target, "Marmota_marmota_marmota.marMar2.1")


@pytest.mark.skipif(not travis, reason="slow")
def test_download_annotation(url):
    out_dir = os.getcwd()
    target = (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/"
        "GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.gff.gz"
    )
    localname = "my_annot"
    with TemporaryDirectory(dir=out_dir) as tmpdir:
        url.download_annotation(
            url="string",
            genomes_dir=tmpdir,
            localname=localname,
            **{"to_annotation": target},
        )

        # check download_and_generate_annotation output
        fname = os.path.join(tmpdir, localname, localname + ".annotation.gtf")
        validate_annot(fname, "gtf")

        fname = os.path.join(tmpdir, localname, localname + ".annotation.bed")
        validate_annot(fname, "bed")

        # check attempt_download_and_report_back output
        readme = os.path.join(tmpdir, localname, "README.txt")
        metadata, lines = genomepy.utils.read_readme(readme)
        assert metadata["annotation url"] == target
