import os
from shutil import copyfile
from tempfile import TemporaryDirectory


def test_ncbiprovider__init__(ncbi):
    assert ncbi.name == "NCBI"
    assert ncbi.taxid_fields == ["species_taxid", "taxid"]


def test__get_genomes(ncbi):
    assert isinstance(ncbi.genomes, dict)
    assert "ASM2732v1" in ncbi.genomes
    genome = ncbi.genomes["ASM2732v1"]
    assert isinstance(genome, dict)
    for field in ncbi.accession_fields + ncbi.taxid_fields + ncbi.description_fields:
        assert field in genome
    assert genome["species_taxid"] == "2097"
    assert genome["taxid"] == "243273"


def test_genome_info_tuple(ncbi):
    t = ncbi._genome_info_tuple("ASM2732v1")
    assert isinstance(t, tuple)
    assert t[2:4] == ("Mycoplasma genitalium G37", "2097")


def test_get_genome_download_link(ncbi):
    link = ncbi.get_genome_download_link("ASM2732v1")
    assert (
        link
        == "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/"
        + "325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz"
    )


def test__post_process_download(ncbi):
    name = "ASM2732v1"
    localname = "tmp"
    out_dir = os.getcwd()
    with TemporaryDirectory(dir=out_dir) as tmpdir:
        # copy fa file for unmasking
        g = os.path.join(tmpdir, localname + ".fa")
        copyfile("tests/data/gap.fa", g)

        ncbi._post_process_download(
            name=name, localname=localname, out_dir=tmpdir, mask="hard"
        )
        assert os.path.exists(g)
        with open(g) as f:
            ln = 0
            for line in f:
                if line.startswith(">"):
                    # test sequence names (should remain identical here)
                    names = line[1:].strip().split(" ")
                    assert names[0] == names[1]
                elif ln == 9:
                    # test masking
                    assert line.strip() == "ANNNNNNNA"
                ln += 1


def test_get_annotation_download_link(ncbi):
    link = ncbi.get_annotation_download_link("ASM2732v1")
    assert (
        link
        == "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/"
        + "000/027/325/GCF_000027325.1_ASM2732v1/"
        + "GCF_000027325.1_ASM2732v1_genomic.gff.gz"
    )
