import os

import pandas as pd
import pytest

import genomepy
import genomepy.annotation.sanitize
import genomepy.files
import genomepy.utils
from genomepy.annotation import query_mygene


def test__get_name_and_dir():
    # working inputs
    t_dicts = [
        {
            "name": "GRCz11",
            "genomes_dir": "tests/data",
            "expected_name": "GRCz11",
            "expected_dir": "tests/data/GRCz11",
        },
        {
            "name": "tests/data/GRCz11",
            "genomes_dir": None,
            "expected_name": "GRCz11",
            "expected_dir": "tests/data/GRCz11",
        },
        {
            "name": "tests/data/GRCz11/GRCz11.annotation.bed",
            "genomes_dir": None,
            "expected_name": "GRCz11",
            "expected_dir": "tests/data/GRCz11",
        },
        # works with warnings
        {
            "name": "empty",
            "genomes_dir": "tests/data",
            "expected_name": "empty",
            "expected_dir": "tests/data/empty",
        },
    ]
    for td in t_dicts:
        name, genome_dir = genomepy.annotation._get_name_and_dir(
            td["name"], td["genomes_dir"]
        )
        assert name == td["expected_name"], td
        assert genome_dir == genomepy.utils.cleanpath(td["expected_dir"]), td

    # does not exist
    with pytest.raises(FileNotFoundError):
        genomepy.annotation._get_name_and_dir("tests/data/GRCz11/what.who")
    # no bed/gtf/fa
    with pytest.raises(NotImplementedError):
        genomepy.annotation._get_name_and_dir("tests/data/GRCz11/README.txt")


def test_annotation_init(caplog, annot):
    assert annot.genome_file.endswith("data/regexp/regexp.fa")
    assert annot.annotation_gtf_file.endswith("data/regexp/regexp.annotation.gtf")
    assert annot.annotation_bed_file.endswith("data/regexp/regexp.annotation.bed")

    assert annot.annotation_contigs == ["chrM"]  # from BED
    assert len(annot.genome_contigs) == 17  # from sizes
    assert isinstance(annot.bed, pd.DataFrame)
    assert isinstance(annot.gtf, pd.DataFrame)
    assert annot.genes("bed") == ["NP_059343.1"]

    # >1 GTF files: take unzipped
    g1 = "tests/data/data.annotation.gtf"
    with genomepy.files._open(g1, "w") as fa:
        fa.write("1")
    g2 = "tests/data/data.annotation.gtf.gz"
    with genomepy.files._open(g2, "w") as fa:
        fa.write("2")
    a = genomepy.Annotation("data", genomes_dir="tests")
    assert a.annotation_gtf_file.endswith(g1)
    genomepy.utils.rm_rf(g1)
    genomepy.utils.rm_rf(g2)

    # not enough GTF files
    genomepy.Annotation("empty", genomes_dir="tests/data")
    assert "Could not find 'empty.annotation.bed" in caplog.text
    assert "Could not find 'empty.annotation.gtf" in caplog.text

    # Genome doesn't exist
    with pytest.raises(FileNotFoundError):
        genomepy.Annotation("never_existed", genomes_dir="tests/data")


def test_custom_annotation():
    for fname in [
        "tests/data/custom.annotation.bed",
        "tests/data/custom.annotation.bed.gz",
    ]:
        a = genomepy.Annotation(name=fname, genomes_dir="tests/data")
        assert a.bed.shape[0] == 10

    for fname in [
        "tests/data/custom.annotation.gtf",
        "tests/data/custom.annotation.gtf.gz",
    ]:
        a = genomepy.Annotation(name=fname, genomes_dir="tests/data")
        assert a.gtf.shape[0] == 45

    with pytest.raises(NotImplementedError):
        a = genomepy.Annotation(name="tests/data/regions.txt", genomes_dir="tests/data")


def test_named_gtf():
    a = genomepy.Annotation("sacCer3", genomes_dir="tests/data")
    df = a.named_gtf
    assert df.index.name == "gene_name"
    assert str(a.gtf.index.dtype) == "int64"
    assert str(df.index.dtype) == "object"
    assert set(df.at["YDL248W", "seqname"]) == {"chrIV"}


def test_genes():
    a = genomepy.Annotation("GRCz11", genomes_dir="tests/data")
    g = a.genes("bed")
    assert isinstance(g, list)
    g = a.genes("gtf")
    assert isinstance(g, list)


def test_gene_coords(caplog):
    a = genomepy.Annotation("sacCer3", genomes_dir="tests/data")

    bed_genes = a.genes("bed")[0:10]
    c = a.gene_coords(bed_genes, "bed")
    assert list(c.shape) == [10, 5]
    assert c.columns.to_list() == ["chrom", "start", "end", "name", "strand"]

    gtf_genes = a.genes("gtf")[0:10]
    c = a.gene_coords(gtf_genes, "gtf")
    assert list(c.shape) == [10, 5]
    assert c.columns.to_list() == ["seqname", "start", "end", "gene_name", "strand"]

    # mismatched gene names!
    _ = a.gene_coords(["what?"], "bed")
    assert "No genes found." in caplog.text

    _ = a.gene_coords(["YDL200C_mRNA", "what?"], "bed")
    assert "Only 50% of genes was found." in caplog.text


def test_map_locations():
    a = genomepy.Annotation("sacCer3", genomes_dir="tests/data")

    # BED + Ensembl
    ens = a.map_locations(annot=a.bed.head(1), to="Ensembl")
    assert ens.chrom.to_list() == ["IV"]

    # GTF + NCBI
    ncbi = a.map_locations(annot=a.gtf.head(1), to="NCBI")
    assert ncbi.seqname.to_list() == ["IV"]

    # custom dataframe (already indexed)
    df = a.bed.head(1).set_index("start")
    cus = a.map_locations(annot=df, to="Ensembl")
    assert ens.chrom.to_list() == ["IV"]
    assert cus.index.name == "start"


def test_filter_regex():
    # match
    gtf_file = "tests/data/regexp/regexp.annotation.gtf"
    with open(gtf_file, "w") as f:
        f.write("chr2\tcool_gene\n")
        f.write("chromosome3\tboring_gene\n")
    a = genomepy.annotation.Annotation("regexp", genomes_dir="tests/data")
    gtf_df = pd.read_csv(gtf_file, sep="\t", header=None)

    df = a.filter_regex(gtf_df, regex=".*2")
    assert isinstance(df, pd.DataFrame)
    assert df[0].to_list() == ["chr2"]

    # invert match
    df = a.filter_regex(gtf_df, regex=".*2", invert_match=True)
    assert isinstance(df, pd.DataFrame)
    assert df[0].to_list() == ["chromosome3"]

    genomepy.utils.rm_rf(gtf_file)


def test_gtf_dict():
    a = genomepy.annotation.Annotation("GRCz11", genomes_dir="tests/data")
    gid2gname = a.gtf_dict("gene_id", "gene_name")
    assert gid2gname["ENSDARG00000103202"] == "CR383668.1"
    gid2gname = a.gtf_dict("gene_id", "gene_name", False)
    assert gid2gname["ENSDARG00000103202"] == ["CR383668.1"]

    seq2tid = a.gtf_dict("seqname", "transcript_id")
    assert seq2tid["1"] == ["ENSDART00000171868"]
    #
    # with pytest.raises(IndexError):
    #     a.gtf_dict("gene_name", "transcript_id", annot="bed")


def test_lengths():
    a = genomepy.annotation.Annotation("GRCz11", genomes_dir="tests/data")

    lengths = a.lengths("gene_name")
    assert str(lengths.dtype) == "uint32"
    expected_genes = a.named_gtf[a.named_gtf.feature == "exon"].index.unique()
    assert sorted(expected_genes) == sorted(lengths.index)
    assert lengths.loc["CR383668.1"] == 420

    lengths = a.lengths("gene_id")
    assert lengths.index[0].startswith("ENSDARG")
    assert lengths.loc["ENSDARG00000103202"] == 420

    lengths = a.lengths("transcript_name")
    assert lengths.loc["CR383668.1-201"] == 420

    lengths = a.lengths("transcript_id")
    assert lengths.index[0].startswith("ENSDART")
    assert lengths.loc["ENSDART00000159919"] == 420


# annotation.utils.py


def test_count_columns():
    cols = genomepy.annotation.utils.count_columns(
        "tests/data/sacCer3/sacCer3.annotation.gtf"
    )
    assert cols == 9

    cols = genomepy.annotation.utils.count_columns(
        "tests/data/sacCer3/sacCer3.annotation.bed"
    )
    assert cols == 12


def test_read_annot():
    df = genomepy.annotation.utils.read_annot(
        "tests/data/sacCer3/sacCer3.annotation.bed"
    )
    assert len(df.columns) == 12
    assert list(df.columns).__eq__(list(genomepy.annotation.utils.BED12_FORMAT))


def test_write_annot():
    for ext in ["bed", "gtf"]:
        df1 = genomepy.annotation.utils.read_annot(
            f"tests/data/sacCer3/sacCer3.annotation.{ext}"
        )
        genomepy.annotation.utils.write_annot(
            df1, f"tests/data/sacCer3/test.annotation.{ext}"
        )
        df2 = genomepy.annotation.utils.read_annot(
            f"tests/data/sacCer3/test.annotation.{ext}"
        )

        os.remove(f"tests/data/sacCer3/test.annotation.{ext}")
        assert df1.equals(df2)


def test_generate_annot():
    bed_file = "tests/data/data.annotation.bed.gz"
    with genomepy.files._open(bed_file, "w") as f:
        f.write("this is the old bed\n")

    gtf_file = "tests/data/data.annotation.gtf.gz"
    with genomepy.files._open(gtf_file, "w") as f:
        f.write(
            "chr1\tgenomepy\texon\t1\t100\t.\t+\t.\t"
            'gene_id "ENSGP1234"; transcript_id "GP_1234.1";  gene_name "GP";\n'
        )
    a = genomepy.annotation.Annotation("data", genomes_dir="tests")

    # BED file before & after
    with genomepy.files._open(a.annotation_bed_file, "r") as f:
        lines = f.readlines()
    assert lines == ["this is the old bed\n"]

    genomepy.annotation.utils.generate_annot(
        a.annotation_gtf_file, a.annotation_bed_file, True
    )
    with genomepy.files._open(a.annotation_bed_file, "r") as f:
        lines = f.readlines()
    assert lines == ["chr1\t0\t100\tGP_1234.1\t0\t+\t100\t100\t0\t1\t100,\t0,\n"]
    assert len(lines) == 1

    # unzipping and rezipping correctly
    assert a.annotation_gtf_file.endswith(".gz")
    assert a.annotation_bed_file.endswith(".gz")

    # does not overwrite by default
    with pytest.raises(FileExistsError):
        genomepy.annotation.utils.generate_annot(
            a.annotation_gtf_file, a.annotation_bed_file
        )

    genomepy.utils.rm_rf(gtf_file)
    genomepy.utils.rm_rf(bed_file)


def test__parse_annot():
    a = genomepy.Annotation("sacCer3", genomes_dir="tests/data")
    df = genomepy.annotation.utils._parse_annot(a, "bed")
    assert df.equals(a.bed)
    df = genomepy.annotation.utils._parse_annot(a, "gtf")
    assert df.equals(a.gtf)
    df = genomepy.annotation.utils._parse_annot(a, a.bed)
    assert df.equals(a.bed)

    with pytest.raises(ValueError):
        genomepy.annotation.utils._parse_annot(a, 1)


# annotation.sanitize.py


def test_match_contigs():
    # nothing to work with
    a = genomepy.Annotation("sacCer3", genomes_dir="tests/data")
    cd = genomepy.annotation.sanitize._match_contigs(a)
    assert cd is None

    # one missing contig, one fixable contig
    a = genomepy.Annotation("sanitize", genomes_dir="tests/data")
    before = a.gtf
    cd = genomepy.annotation.sanitize._match_contigs(a)
    after = a.gtf

    assert cd == {"NC_007112.7": "1"}
    assert a.genome_contigs == ["1"]
    assert list(before.seqname.unique()) == ["NC_007112.7", "NC_002333.2"]
    assert list(after.seqname.unique()) == ["1", "NC_002333.2"]


def test_filter_contigs():
    a = genomepy.Annotation("sacCer3", genomes_dir="tests/data")
    assert "chrV" not in a.bed.chrom.unique()

    # add a chromosome to the BED not present in the genome
    a.bed.at[0, "chrom"] = "chrV"
    assert "chrV" in a.bed.chrom.unique()

    missing_contigs = genomepy.annotation.sanitize._filter_contigs(a)
    assert missing_contigs == {"chrV"}
    assert "chrV" not in a.bed.chrom.unique()


def test_document_sanitizing():
    class HasReadme:
        readme_file = "tests/data/test.log"

        def __init__(self):
            with open(self.readme_file, "w"):
                pass

        def _check_property(self):
            pass

    a = HasReadme()
    cd = {"1": "a", "2": "b"}
    mc = ["1", "2", "3", "4", "5"]
    genomepy.annotation.sanitize._document_sanitizing(a, cd, mc)
    md, lines = genomepy.files.read_readme(a.readme_file)
    genomepy.utils.rm_rf(a.readme_file)

    status = md["sanitized annotation"]
    assert status == "2 contigs were renamed. 5 contigs were removed (see below)."
    assert lines == [
        "",
        "The following contigs were filtered out of the gene annotation:",
        "1, 2, 3, 4, 5",
    ]


def test_sanitize():
    a = genomepy.Annotation("sanitize", genomes_dir="tests/data")

    a.sanitize(overwrite=False)
    assert a.gtf.shape == (4, 9)
    assert list(set(a.gtf.seqname)) == ["1"]


# annotation.mygene.py


def test_map_genes():
    a = genomepy.Annotation("GRCz11", genomes_dir="tests/data")

    bed = a.bed.head()
    transcript_ids = bed.name.to_list()
    assert transcript_ids[0] == "ENSDART00000159919"

    # transcript to gene
    res = a.map_genes(field="ensembl.gene", annot=bed)
    genes = res.name.to_list()
    assert genes[0] == "ENSDARG00000103202"

    # # transcript to symbol
    # res = a.map_genes(field="symbol", df=bed)
    # symbol = res.name.to_list()
    # assert symbol[0] == "CR383668.1"

    # refseq hits & subtypes
    protein = a.map_genes(field="refseq", product="protein", annot=bed)
    assert protein.name.to_list()[0].startswith("NP_")
    # rna = a.map_genes(field="refseq", product="rna", df=bed)
    # assert rna.name.to_list()[0].startswith("NM_")
    # assert rna.shape == protein.shape


def test_query_mygene():
    result = query_mygene(["ENST00000449992"], 9606, "symbol")
    assert result.loc["ENST00000449992", "symbol"] == "TP63"


def test_parse_mygene_input():
    # wrong product
    with pytest.raises(ValueError):
        genomepy.annotation.mygene._parse_mygene_input("", "illegal!")

    # wrong field
    with pytest.raises(ValueError):
        genomepy.annotation.mygene._parse_mygene_input("illegal!", None)

    # correct input
    gf, p = genomepy.annotation.mygene._parse_mygene_input(
        "refseq.translation.rna", None
    )
    assert gf == "refseq.translation.rna" and p is None
    gf, p = genomepy.annotation.mygene._parse_mygene_input("NAME", "RNA")
    assert gf == "name" and p == "rna"


# def test_ensembl_genome_info():
#     # Works
#     a = genomepy.Annotation("sacCer3", genomes_dir="tests/data")
#     egi = genomepy.annotation.mygene.ensembl_genome_info(a)
#     assert egi == ("R64-1-1", "GCA_000146045.2", 4932)
#     # genomepy.utils.rm_rf(os.path.join(a.genome_dir, "assembly_report.txt"))
#
#     # No readme
#     a = genomepy.Annotation("regexp", genomes_dir="tests/data")
#     readme_file = os.path.join(a.genome_dir, "README.txt")
#     with pytest.raises(FileNotFoundError):
#         genomepy.annotation.mygene.ensembl_genome_info(a)
#
#     # Empty readme
#     with open(readme_file, "w") as f:
#         f.write("\n")
#     a = genomepy.Annotation("regexp", genomes_dir="tests/data")
#     egi = genomepy.annotation.mygene.ensembl_genome_info(a)
#     assert egi is None
#     genomepy.utils.rm_rf(readme_file)
