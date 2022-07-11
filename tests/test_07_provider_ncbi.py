import os
from shutil import copyfile
from tempfile import TemporaryDirectory

import pandas as pd

import genomepy


def test_ncbiprovider(ncbi):
    assert ncbi.name == "NCBI"
    assert ncbi.taxid_fields == ["species_taxid", "taxid"]


def test_genome_info_tuple(ncbi):
    t = ncbi._genome_info_tuple("ASM2732v1", size=True)
    assert isinstance(t, tuple)
    assert t[:-1] == (
        "ASM2732v1",
        "GCF_000027325.1",
        2097,
        True,
        "Mycoplasma genitalium G37",
        580076,
    )


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

        ncbi._post_process_download(name=name, fname=g, out_dir=tmpdir, mask="hard")
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

        asm_report = os.path.join(tmpdir, "assembly_report.txt")
        assert os.path.exists(asm_report)


def test_get_annotation_download_link(ncbi):
    links = ncbi.get_annotation_download_links("ASM2732v1")
    assert links[0].endswith("GCF_000027325.1_ASM2732v1_genomic.gff.gz")


def test__get_genomes(ncbi):
    assert isinstance(ncbi.genomes, dict)
    assert "ASM2732v1" in ncbi.genomes
    genome = ncbi.genomes["ASM2732v1"]
    assert isinstance(genome, dict)
    for field in ncbi.accession_fields + ncbi.taxid_fields + ncbi.description_fields:
        assert field in genome
    assert genome["species_taxid"] == "2097"
    assert genome["taxid"] == "243273"


def test_download_assembly_report():
    assembly_report = "tests/data/sacCer3/assembly_report.txt"
    genomepy.providers.download_assembly_report("GCA_000146045", assembly_report)
    report = pd.read_csv(assembly_report, sep="\t", comment="#")

    assert isinstance(report, pd.DataFrame)
    assert list(report.columns) == genomepy.providers.ncbi.ASM_FORMAT
