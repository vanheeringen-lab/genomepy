import genomepy
import gzip
import os
import pytest
import shutil

from tempfile import TemporaryDirectory


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
    return genomepy.provider.NcbiProvider()


def test_ncbiprovider__init__(p):
    p2 = genomepy.provider.ProviderBase().create("NCBI")
    assert p2.name == p.name == "NCBI"
    assert p.taxid_fields == ["species_taxid", "taxid"]


def test__get_genomes(p):
    assert isinstance(p.genomes, dict)
    assert "ASM2732v1" in p.genomes
    genome = p.genomes["ASM2732v1"]
    assert isinstance(genome, dict)
    for field in p.accession_fields + p.taxid_fields + p.description_fields:
        assert field in genome
    assert genome["species_taxid"] == "2097"
    assert genome["taxid"] == "243273"


def test_genome_info_tuple(p):
    t = p._genome_info_tuple("ASM2732v1")
    assert isinstance(t, tuple)
    assert t[2:4] == ("Mycoplasma genitalium G37", "2097")


def test_get_genome_download_link(p):
    link = p.get_genome_download_link("ASM2732v1")
    assert (
        link
        == "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/"
        + "325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz"
    )


def test__post_process_download(p):
    name = "ASM2732v1"
    localname = "tmp"
    out_dir = os.getcwd()
    with TemporaryDirectory(dir=out_dir) as tmpdir:
        # copy fa file for unmasking
        g = os.path.join(tmpdir, localname + ".fa")
        shutil.copyfile("tests/data/gap.fa", g)

        p._post_process_download(
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


def test_get_annotation_download_link(p):
    link = p.get_annotation_download_link("ASM2732v1")
    assert (
        link
        == "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/"
        + "000/027/325/GCF_000027325.1_ASM2732v1/"
        + "GCF_000027325.1_ASM2732v1_genomic.gff.gz"
    )
