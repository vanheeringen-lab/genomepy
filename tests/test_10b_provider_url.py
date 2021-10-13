import pytest

import genomepy.providers.url
import genomepy.utils


def test_urlprovider(url):
    assert url.name == "URL"
    assert url.genomes == {}


def test_genome_taxid(url):
    assert url.genome_taxid(None) is None


def test_assembly_accession(url):
    assert url.assembly_accession(None) is None


def test_search(url):
    result = url.search("???")
    with pytest.raises(StopIteration):
        assert next(result)


def test__genome_info_tuple(url):
    assert url._genome_info_tuple(None) is tuple()


def test__check_name(url):
    assert url._check_name(None) is None


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


def test_fuzzy_annotation_search():
    search_name = "mygenome"
    search_list = [
        "nothing.gtf",
        "mygenome.genes.gff3",
        "asdasdasdmygenomeasdasdasd.gtf",
        "mygenome.zipped.gtf.gz",
    ]
    expected = [
        "mygenomeasdasdasd.gtf",
        "mygenome.zipped.gtf.gz",
        "mygenome.genes.gff3",
    ]
    # f"{search_name}.*?\.{ext}3?(\.gz)?"
    hits = genomepy.providers.url.fuzzy_annotation_search(search_name, search_list)
    assert hits == expected


def test_search_url_for_annotations():
    target = "http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_genome.fa.gz"
    links = genomepy.providers.url.search_url_for_annotations(target)
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
    links = genomepy.providers.url.search_url_for_annotations(target)
    expected = [
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/"
        + "GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.gtf.gz",
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/"
        + "GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.gff.gz",
    ]
    assert links == expected

    # no annot file
    target = (
        "http://ftp.ensembl.org/pub/release-100/fasta/marmota_marmota_marmota/"
        "dna/Marmota_marmota_marmota.marMar2.1.dna.toplevel.fa.gz"
    )
    links = genomepy.providers.url.search_url_for_annotations(target)
    assert links == []
