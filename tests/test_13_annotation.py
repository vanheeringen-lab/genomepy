import os

import pandas as pd
import pytest

import genomepy
import genomepy.files
import genomepy.utils
from genomepy.annotation import query_mygene


def test_annotation_init(caplog, annot):
    assert annot.genome_file.endswith("data/regexp/regexp.fa")
    assert annot.annotation_gtf_file.endswith("data/regexp/regexp.annotation.gtf")
    assert annot.annotation_bed_file.endswith("data/regexp/regexp.annotation.bed")

    assert annot.annotation_contigs == ["chrM"]  # from BED
    assert len(annot.genome_contigs) == 17  # from sizes
    assert isinstance(annot.bed, pd.DataFrame)
    assert isinstance(annot.gtf, pd.DataFrame)
    assert annot.genes() == ["NP_059343.1"]

    # >1 GTF files: take unzipped
    g1 = "tests/data/data.annotation.gtf"
    with genomepy.files._open(g1, "w") as fa:
        fa.write("1")
    g2 = "tests/data/data.annotation.gtf.gz"
    with genomepy.files._open(g2, "w") as fa:
        fa.write("2")
    a = genomepy.Annotation(name="data", genomes_dir="tests")
    assert a.annotation_gtf_file.endswith(g1)
    genomepy.utils.rm_rf(g1)
    genomepy.utils.rm_rf(g2)

    # not enough GTF files
    genomepy.Annotation(name="never_existed", genomes_dir="tests/data")
    assert "Could not find 'never_existed.annotation.bed" in caplog.text
    assert "Could not find 'never_existed.annotation.gtf" in caplog.text


def test_named_gtf():
    a = genomepy.Annotation("sacCer3", "tests/data")
    df = a.named_gtf
    assert df.index.name == "gene_name"
    assert str(a.gtf.index.dtype) == "int64"
    assert str(df.index.dtype) == "object"
    assert set(df.at["YDL248W", "seqname"]) == {"chrIV"}


def test_genes(annot):
    g = annot.genes()
    assert isinstance(g, list)
    g = annot.genes("gtf")
    assert isinstance(g, list)


def test_gene_coords(caplog):
    a = genomepy.Annotation("sacCer3", "tests/data")

    bed_genes = a.genes()[0:10]
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
    a = genomepy.Annotation("sacCer3", "tests/data")

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
    a = genomepy.annotation.Annotation("regexp", "tests/data")
    gtf_df = pd.read_csv(gtf_file, sep="\t", header=None)

    df = a.filter_regex(gtf_df, regex=".*2")
    assert isinstance(df, pd.DataFrame)
    assert df[0].to_list() == ["chr2"]

    # invert match
    df = a.filter_regex(gtf_df, regex=".*2", invert_match=True)
    assert isinstance(df, pd.DataFrame)
    assert df[0].to_list() == ["chromosome3"]

    genomepy.utils.rm_rf(gtf_file)


# annotation.utils.py


def test_count_columns():
    pass  # TODO


def test_read_annot():
    pass  # TODO


def test_write_annot():
    pass  # TODO


def test_generate_annot():
    # bed_file = "tests/data/data.annotation.bed.gz"
    # with genomepy.annotation._open(bed_file, "w") as f:
    #     f.write("this is the old bed\n")
    #
    # gtf_file = "tests/data/data.annotation.gtf.gz"
    # with genomepy.annotation._open(gtf_file, "w") as f:
    #     f.write(
    #         "chr1\tgenomepy\texon\t1\t100\t.\t+\t.\t"
    #         'gene_id "ENSGP1234"; transcript_id "GP_1234.1";  gene_name "GP";\n'
    #     )
    # a = genomepy.annotation.Annotation("data", "tests")
    #
    # # BED file before & after
    # with genomepy.annotation._open(a.annotation_bed_file, "r") as f:
    #     lines = f.readlines()
    # assert lines == ["this is the old bed\n"]
    # a.bed_from_gtf()
    # with genomepy.annotation._open(a.annotation_bed_file, "r") as f:
    #     lines = f.readlines()
    # assert lines == ["chr1\t0\t100\tGP_1234.1\t0\t+\t100\t100\t0\t1\t100,\t0,\n"]
    # assert len(lines) == 1
    #
    # # unzipping and rezipping correctly
    # assert a.annotation_gtf_file.endswith(".gz")
    # assert a.annotation_bed_file.endswith(".gz")
    #
    # genomepy.utils.rm_rf(gtf_file)
    # genomepy.utils.rm_rf(bed_file)
    pass  # TODO


def test__parse_annot():
    # a = genomepy.Annotation("sacCer3", "tests/data")
    # df = a._parse_annot("bed")
    # assert df.equals(a.bed)
    # df = a._parse_annot("gtf")
    # assert df.equals(a.gtf)
    # df = a._parse_annot(a.bed)
    # assert df.equals(a.bed)
    # assert a._parse_annot(1) is None
    pass  # TODO


# annotation.sanitize.py


def test_match_contigs():
    # a = genomepy.Annotation("sacCer3", "tests/data")
    # cd = genomepy.annotation.sanitize.match_contigs(a)
    pass  # TODO


def test_filter_contigs():
    pass  # TODO


def test_document_sanitizing():
    pass  # TODO


def test_sanitize():
    pass  # TODO


# def test_filter_genome_contigs():
#     # chr2 found, chromosome3 missing from genome
#     gtf_file = "tests/data/regexp/regexp.annotation.gtf"
#     with open(gtf_file, "w") as f:
#         f.write("chr2\tcool_gene\n")
#         f.write("chromosome3\tboring_gene\n")
#     a = genomepy.annotation.Annotation("regexp", "tests/data")
#     genomepy.files.generate_fa_sizes(a.genome_file, a.genome_file + ".sizes")
#     a.sizes_file = a.genome_file + ".sizes"
#
#     missing_contigs = a._filter_genome_contigs()
#     assert missing_contigs == ["chromosome3"]
#
#     genomepy.utils.rm_rf(gtf_file)
#     genomepy.utils.rm_rf(a.genome_file + ".sizes")
#
#
# def test__conforming_index():
#     # genome header found in annotation
#     bed_file = "tests/data/regexp/regexp.annotation.bed"
#     with open(bed_file, "w") as f:
#         f.write(
#             """chr1\t15307\t16448\tNP_059343.1\t42\t+\t15307\t16448\t0\t1\t1141,\t0,"""
#         )
#
#     a = genomepy.annotation.Annotation("regexp", "tests/data")
#     i = a._conforming_index()
#     assert i == 0
#
#     # genome header not found in annotation
#     bed_file = "tests/data/regexp/regexp.annotation.bed"
#     with open(bed_file, "w") as f:
#         f.write(
#             """chrM\t15307\t16448\tNP_059343.1\t42\t+\t15307\t16448\t0\t1\t1141,\t0,"""
#         )
#
#     a = genomepy.annotation.Annotation("regexp", "tests/data")
#     i = a._conforming_index()
#     assert i == -1
#
#     genomepy.utils.rm_rf(bed_file)
#
#
# def test__contig_conversion_dict():
#     gtf_file = "tests/data/data.annotation.gtf"
#     with open(gtf_file, "w") as f:
#         f.write("\n")
#     bed_file = "tests/data/data.annotation.bed"
#     with open(bed_file, "w") as f:
#         f.write("\n")
#     genome_file = "tests/data/data.fa"
#     with open(genome_file, "w") as f:
#         f.write(">this matches chr2\n" ">that matches chr3\n")
#     a = genomepy.annotation.Annotation("data", "tests")
#
#     conversion_dict, duplicate_contigs = a._contig_conversion_dict(2)
#     assert not duplicate_contigs
#     assert conversion_dict["chr2"] == "this"
#     assert conversion_dict["chr3"] == "that"
#
#     # assume "matches" is the contig name
#     ids, duplicate_contigs = a._contig_conversion_dict(1)
#     assert "matches" in duplicate_contigs
#     assert ids["matches"] == "that"
#
#     genomepy.utils.rm_rf(gtf_file)
#     genomepy.utils.rm_rf(bed_file)
#     genomepy.utils.rm_rf(genome_file)
#
#
# def test__conform_gtf():
#     bed_file = "tests/data/data.annotation.bed"
#     with open(bed_file, "w") as f:
#         f.write("\n")
#
#     # filter missing contigs
#     conversion_dict = {"chr1": "chromosome_one"}
#     gtf_file = "tests/data/data.annotation.gtf"
#     with open(gtf_file, "w") as f:
#         f.write("chr1\t1\t2\n")
#         f.write("chr2\t1\t2\n")
#     a = genomepy.annotation.Annotation("data", "tests")
#     missing_contigs = a._conform_gtf(conversion_dict)
#     assert missing_contigs == ["chr2"]
#     with open(gtf_file) as f:
#         lines = f.readlines()
#     assert lines == ["chromosome_one\t1\t2\n"]
#
#     # keep missing contigs
#     with open(gtf_file, "w") as f:
#         f.write("chr1\t1\t2\n")
#         f.write("chr2\t1\t2\n")
#     a = genomepy.annotation.Annotation("data", "tests")
#     missing_contigs = a._conform_gtf(conversion_dict, filter_contigs=False)
#     assert missing_contigs == ["chr2"]
#     with open(gtf_file) as f:
#         lines = f.readlines()
#     assert lines == ["chromosome_one\t1\t2\n", "chr2\t1\t2\n"]
#
#     genomepy.utils.rm_rf(gtf_file)
#     genomepy.utils.rm_rf(bed_file)
#
#
# def test_sanitize(caplog):
#     tmp_dir = mkdtemp(dir="tests")
#     genome = "testgenome"
#     genomepy.utils.mkdir_p(os.path.join(tmp_dir, genome))
#
#     bed_file = os.path.join(tmp_dir, genome, genome + ".annotation.bed")
#     gtf_file = os.path.join(tmp_dir, genome, genome + ".annotation.gtf")
#     genome_file = os.path.join(tmp_dir, genome, genome + ".fa")
#     sizes_file = os.path.join(tmp_dir, genome, genome + ".fa.sizes")
#     readme_file = os.path.join(tmp_dir, genome, "README.txt")
#
#     # no genome
#     a = genomepy.annotation.Annotation(genome, tmp_dir)
#     a.sanitize()
#     assert "A genome is required for sanitizing!" in caplog.text
#
#     # conforming, filtering off
#     with open(bed_file, "w") as f:
#         f.write("chr1\t0\t100\n" "chr2\t0\t100\n")
#     with open(gtf_file, "w") as f:
#         f.write(
#             "chr1\tgenomepy\texon\t1\t100\t.\t+\t.\t"
#             'gene_id "ENSGP1234"; transcript_id "GP_1234.1";  gene_name "GP1";\n'
#             "chr2\tgenomepy\texon\t1\t100\t.\t+\t.\t"
#             'gene_id "ENSGP1235"; transcript_id "GP_1235.1";  gene_name "GP2";\n'
#         )
#     with open(genome_file, "w") as f:
#         f.write(">chr1\n")
#         f.write("ATCGATCG\n")
#     with open(readme_file, "w") as f:
#         f.write("\n")
#     genomepy.files.generate_fa_sizes(genome_file, sizes_file)
#     a = genomepy.annotation.Annotation(genome, tmp_dir)
#     a.sanitize(match_contigs=False, filter_contigs=False)
#
#     # expect no changes
#     with open(gtf_file) as f:
#         lines = f.readlines()
#     assert lines == [
#         'chr1\tgenomepy\texon\t1\t100\t.\t+\t.\tgene_id "ENSGP1234"; transcript_id "GP_1234.1";  gene_name "GP1";\n',
#         'chr2\tgenomepy\texon\t1\t100\t.\t+\t.\tgene_id "ENSGP1235"; transcript_id "GP_1235.1";  gene_name "GP2";\n',
#     ]
#     with open(bed_file) as f:
#         lines = f.readlines()
#     assert lines == ["chr1\t0\t100\n", "chr2\t0\t100\n"]
#     metadata, _ = genomepy.files.read_readme(a.readme_file)
#     assert metadata["sanitized annotation"] == "contigs match but not filtered"
#
#     # conforming, filtering on
#     a.sanitize(match_contigs=False, filter_contigs=True)
#
#     # expect chr2 to have been removed
#     with open(gtf_file) as f:
#         lines = f.readlines()
#     assert lines == [
#         'chr1\tgenomepy\texon\t1\t100\t.\t+\t.\tgene_id "ENSGP1234"; transcript_id "GP_1234.1";  gene_name "GP1";\n',
#     ]
#     with open(bed_file) as f:
#         lines = f.readlines()
#     assert lines == [
#         "chr1\t0\t100\tGP_1234.1\t0\t+\t100\t100\t0\t1\t100,\t0,\n",
#     ]
#     metadata, _ = genomepy.files.read_readme(a.readme_file)
#     assert metadata["sanitized annotation"] == "contigs match and filtered"
#
#     # not conforming, no fix possible
#     genomepy.utils.rm_rf(os.path.join(tmp_dir, genome))
#     genomepy.utils.mkdir_p(os.path.join(tmp_dir, genome))
#     with open(bed_file, "w") as f:
#         f.write("chr1\t0\t100\n")
#     with open(gtf_file, "w") as f:
#         f.write("\n")
#     with open(genome_file, "w") as f:
#         f.write(">not_chr1\n")
#         f.write("ATCGATCG\n")
#     with open(readme_file, "w") as f:
#         f.write("\n")
#     genomepy.files.generate_fa_sizes(genome_file, sizes_file)
#     a = genomepy.annotation.Annotation(genome, tmp_dir)
#     a.sanitize()
#     metadata, _ = genomepy.files.read_readme(a.readme_file)
#     assert metadata["sanitized annotation"] == "not possible"
#
#     # not conforming, filtering off
#     genomepy.utils.rm_rf(os.path.join(tmp_dir, genome))
#     genomepy.utils.mkdir_p(os.path.join(tmp_dir, genome))
#     with open(bed_file, "w") as f:
#         f.write("chr1\t0\t100\n" "chr2\t0\t100\n")
#     with open(gtf_file, "w") as f:
#         f.write(
#             "chr1\tgenomepy\texon\t1\t100\t.\t+\t.\t"
#             'gene_id "ENSGP1234"; transcript_id "GP_1234.1";  gene_name "GP1";\n'
#             "chr2\tgenomepy\texon\t1\t100\t.\t+\t.\t"
#             'gene_id "ENSGP1235"; transcript_id "GP_1235.1";  gene_name "GP2";\n'
#         )
#     with open(genome_file, "w") as f:
#         f.write(">this matches chr1\n")
#         f.write("ATCGATCG\n")
#         f.write(">this2 matches chr1\n")
#         f.write("ATCGATCG\n")
#     with open(readme_file, "w") as f:
#         f.write("\n")
#     genomepy.files.generate_fa_sizes(genome_file, sizes_file)
#     a = genomepy.annotation.Annotation(genome, tmp_dir)
#     a.sanitize(match_contigs=True, filter_contigs=False)
#     assert "The genome contains duplicate contig names" in caplog.text
#     with open(bed_file) as f:
#         lines = f.readlines()
#     assert lines == [
#         "this2\t0\t100\tGP_1234.1\t0\t+\t100\t100\t0\t1\t100,\t0,\n",
#         "chr2\t0\t100\tGP_1235.1\t0\t+\t100\t100\t0\t1\t100,\t0,\n",
#     ]
#     metadata, _ = genomepy.files.read_readme(a.readme_file)
#     assert metadata["sanitized annotation"] == "contigs fixed but not filtered"
#
#     # not conforming, filtering off
#     genomepy.utils.rm_rf(os.path.join(tmp_dir, genome))
#     genomepy.utils.mkdir_p(os.path.join(tmp_dir, genome))
#     with open(bed_file, "w") as f:
#         f.write("chr1\t0\t100\n" "chr2\t0\t100\n")
#     with open(gtf_file, "w") as f:
#         f.write(
#             "chr1\tgenomepy\texon\t1\t100\t.\t+\t.\t"
#             'gene_id "ENSGP1234"; transcript_id "GP_1234.1";  gene_name "GP1";\n'
#             "chr2\tgenomepy\texon\t1\t100\t.\t+\t.\t"
#             'gene_id "ENSGP1235"; transcript_id "GP_1235.1";  gene_name "GP2";\n'
#         )
#     with open(genome_file, "w") as f:
#         f.write(">this matches chr1\n")
#         f.write("ATCGATCG\n")
#         f.write(">this2 matches chr1\n")
#         f.write("ATCGATCG\n")
#     with open(readme_file, "w") as f:
#         f.write("\n")
#     genomepy.files.generate_fa_sizes(genome_file, sizes_file)
#     a = genomepy.annotation.Annotation(genome, tmp_dir)
#     a.sanitize(match_contigs=True, filter_contigs=True)
#     assert "The genome contains duplicate contig names" in caplog.text
#     with open(bed_file) as f:
#         lines = f.readlines()
#     assert lines == ["this2\t0\t100\tGP_1234.1\t0\t+\t100\t100\t0\t1\t100,\t0,\n"]
#     metadata, _ = genomepy.files.read_readme(a.readme_file)
#     assert metadata["sanitized annotation"] == "contigs fixed and filtered"
#
#     genomepy.utils.rm_rf(tmp_dir)

# annotation.mygene.py


def test_map_genes():
    a = genomepy.Annotation("GRCz11", "tests/data")

    bed = a.bed.head()
    transcript_ids = bed.name.to_list()
    assert transcript_ids[0] == "ENSDART00000159919"

    # transcript to gene
    res = a.map_genes(gene_field="ensembl.gene", annot=bed)
    genes = res.name.to_list()
    assert genes[0] == "ENSDARG00000103202"

    # # transcript to symbol
    # res = a.map_genes(gene_field="symbol", df=bed)
    # symbol = res.name.to_list()
    # assert symbol[0] == "CR383668.1"

    # refseq hits & subtypes
    protein = a.map_genes(gene_field="refseq", product="protein", annot=bed)
    assert protein.name.to_list()[0].startswith("NP_")
    # rna = a.map_genes(gene_field="refseq", product="rna", df=bed)
    # assert rna.name.to_list()[0].startswith("NM_")
    # assert rna.shape == protein.shape


def test_query_mygene():
    result = query_mygene(["ENST00000449992"], 9606, "symbol")
    assert result.loc["ENST00000449992", "symbol"] == "TP63"


def test_parse_mygene_input():
    pass  # TODO


def test_ensembl_genome_info():
    # Works
    a = genomepy.Annotation("sacCer3", "tests/data")
    egi = genomepy.annotation.mygene.ensembl_genome_info(a)
    assert egi == ("R64-1-1", "GCA_000146045.2", 4932)
    # genomepy.utils.rm_rf(os.path.join(a.genome_dir, "assembly_report.txt"))

    # No readme
    a = genomepy.Annotation("regexp", "tests/data")
    readme_file = os.path.join(a.genome_dir, "README.txt")
    with pytest.raises(FileNotFoundError):
        genomepy.annotation.mygene.ensembl_genome_info(a)

    # Empty readme
    with open(readme_file, "w") as f:
        f.write("\n")
    a = genomepy.Annotation("regexp", "tests/data")
    egi = genomepy.annotation.mygene.ensembl_genome_info(a)
    assert egi is None
    genomepy.utils.rm_rf(readme_file)
