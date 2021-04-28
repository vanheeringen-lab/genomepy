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
    return genomepy.provider.UrlProvider()


def test_ncbiprovider__init__(p):
    p2 = genomepy.provider.ProviderBase().create("URL")
    assert p2.name == p.name == "URL"
    assert p.genomes == {}


def test_genome_taxid(p):
    assert p.genome_taxid({}) == "na"


def test_assembly_accession(p):
    assert p.assembly_accession({}) == "na"


def test_search(p):
    result = p.search("???")
    expected = genomepy.provider.EnsemblProvider().search("???")
    assert isinstance(result, type(expected))
    with pytest.raises(StopIteration):
        assert next(result)


def test_get_genome_download_link(p):
    link = p.get_genome_download_link("url")
    assert link == "url"


def test_get_annotation_download_link(p):
    url = "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_Xenbase.gtf"
    link = p.get_annotation_download_link(None, **{"to_annotation": url})
    assert link == url

    with pytest.raises(TypeError):
        bad_url = "bad_url"
        p.get_annotation_download_link(None, **{"to_annotation": bad_url})

    with pytest.raises(TypeError):
        bad_url = "http://good_url.bad_ext"
        p.get_annotation_download_link(None, **{"to_annotation": bad_url})


def test_search_url_for_annotations(p):
    url = "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_genome.fa.gz"
    links = p.search_url_for_annotations(url, "XENTR_9.1")
    expected = [
        "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_Xenbase.gtf",
        "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_GCA.gff3",
        "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_GCF.gff3",
        "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_Xenbase.gff3",
        "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_Xenbase_longest.gff3",
    ]
    assert len(links) == 5
    assert links == expected

    url = (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/"
        "GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz"
    )
    links = p.search_url_for_annotations(url, "GCF_000027325.1_ASM2732v1")
    expected = [
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/"
        + "GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.gtf.gz",
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/"
        + "GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.gff.gz",
    ]
    assert len(links) == 2
    assert links == expected

    # no annot file
    with pytest.raises(FileNotFoundError):
        url = (
            "http://ftp.ensembl.org/pub/release-100/fasta/marmota_marmota_marmota/"
            "dna/Marmota_marmota_marmota.marMar2.1.dna.toplevel.fa.gz"
        )
        p.search_url_for_annotations(url, "Marmota_marmota_marmota.marMar2.1")


@pytest.mark.skipif(not travis, reason="slow")
def test_download_annotation(p):
    out_dir = os.getcwd()
    annot_url = (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/"
        "GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.gff.gz"
    )
    localname = "my_annot"
    with TemporaryDirectory(dir=out_dir) as tmpdir:
        p.download_annotation(
            url="string",
            genomes_dir=tmpdir,
            localname=localname,
            **{"to_annotation": annot_url},
        )

        # check download_and_generate_annotation output
        fname = os.path.join(tmpdir, localname, localname + ".annotation.gtf.gz")
        validate_gzipped_gtf(fname)

        fname = os.path.join(tmpdir, localname, localname + ".annotation.bed.gz")
        validate_gzipped_bed(fname)

        # check attempt_download_and_report_back output
        readme = os.path.join(tmpdir, localname, "README.txt")
        metadata, lines = genomepy.utils.read_readme(readme)
        assert metadata["annotation url"] == annot_url
