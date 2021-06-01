import os
from tempfile import mkdtemp

import pandas as pd
import pytest

import genomepy
import genomepy.files
import genomepy.utils
from genomepy.annotation import query_mygene


def test_get_column():
    for sep in [",", "\t"]:
        sep_file = "tests/data/test.sep"
        with genomepy.annotation._open(sep_file, "w") as f:
            f.write(f"#comment1{sep}comment2\n")
            f.write(f"c1.1{sep}c2.1\n")
            f.write(f"c1.2{sep}c2.2\n")

        col = genomepy.annotation.get_column(sep_file, sep=sep)
        assert col == ["c1.1", "c1.2"]
        col = genomepy.annotation.get_column(sep_file, n=1, sep=sep, comment_char=None)
        assert col == ["comment2", "c2.1", "c2.2"]
    genomepy.utils.rm_rf(sep_file)


def test_annotation_init(annot):
    assert annot.genome_file.endswith("data/regexp/regexp.fa")
    assert annot.annotation_gtf_file.endswith("data/regexp/regexp.annotation.gtf")
    assert annot.annotation_bed_file.endswith("data/regexp/regexp.annotation.bed")

    assert annot.annotation_contigs == ["chrM"]  # from BED
    assert len(annot.genome_contigs) == 17  # from sizes
    assert isinstance(annot.bed, pd.DataFrame)
    assert isinstance(annot.gtf, pd.DataFrame)
    assert annot.genes == ["NP_059343.1"]

    # too many GTF files: cant choose
    g1 = "tests/data/data.annotation.gtf"
    with genomepy.annotation._open(g1, "w") as fa:
        fa.write("1")
    g2 = "tests/data/data.annotation.gtf.gz"
    with genomepy.annotation._open(g2, "w") as fa:
        fa.write("2")
    with pytest.raises(IndexError):
        genomepy.Annotation(name="data", genomes_dir="tests")
    genomepy.utils.rm_rf(g1)
    genomepy.utils.rm_rf(g2)

    # not enough GTF files
    with pytest.raises(FileNotFoundError):
        genomepy.Annotation(name="never_existed", genomes_dir="tests/data")


def test_gene_coords(annot):
    genes = ["NP_059343.1"]
    coords = annot.gene_coords(genes)
    assert coords.shape[0] == 1  # 1 row per gene


def test_filter_regex():
    bed_file = "tests/data/regexp/regexp.annotation.bed"
    with open(bed_file, "w") as f:
        f.write("\n")

    # match
    gtf_file = "tests/data/regexp/regexp.annotation.gtf"
    with open(gtf_file, "w") as f:
        f.write("chr2\tcool_gene\n")
        f.write("chromosome3\tboring_gene\n")
    a = genomepy.annotation.Annotation("regexp", "tests/data")
    missing_contigs = a.filter_regex(a.annotation_gtf_file, regex=".*2")
    assert missing_contigs == ["chromosome3"]

    # invert match
    gtf_file = "tests/data/regexp/regexp.annotation.gtf"
    with open(gtf_file, "w") as f:
        f.write("chr2\tcool_gene\n")
        f.write("chromosome3\tboring_gene\n")
    missing_contigs = a.filter_regex(
        a.annotation_gtf_file, regex=".*3", invert_match=True
    )
    assert missing_contigs == ["chromosome3"]

    genomepy.utils.rm_rf(gtf_file)
    genomepy.utils.rm_rf(bed_file)


def test_filter_genome_contigs():
    bed_file = "tests/data/regexp/regexp.annotation.bed"
    with open(bed_file, "w") as f:
        f.write("\n")

    # chr2 found, chromosome3 missing from genome
    gtf_file = "tests/data/regexp/regexp.annotation.gtf"
    with open(gtf_file, "w") as f:
        f.write("chr2\tcool_gene\n")
        f.write("chromosome3\tboring_gene\n")
    a = genomepy.annotation.Annotation("regexp", "tests/data")
    genomepy.files.generate_fa_sizes(a.genome_file, a.genome_file + ".sizes")
    missing_contigs = a._filter_genome_contigs()
    assert missing_contigs == ["chromosome3"]

    genomepy.utils.rm_rf(gtf_file)
    genomepy.utils.rm_rf(bed_file)
    genomepy.utils.rm_rf(a.genome_file + ".sizes")


def test__conforming_index():
    # genome header found in annotation
    gtf_file = "tests/data/regexp/regexp.annotation.gtf"
    with open(gtf_file, "w") as f:
        f.write("\n")

    bed_file = "tests/data/regexp/regexp.annotation.bed"
    with open(bed_file, "w") as f:
        f.write("chr1\tcool_gene\n")

    a = genomepy.annotation.Annotation("regexp", "tests/data")
    i = a._conforming_index()
    assert i == 0

    # genome header not found in annotation
    bed_file = "tests/data/regexp/regexp.annotation.bed"
    with open(bed_file, "w") as f:
        f.write("chr2\tcool_gene\n")

    a = genomepy.annotation.Annotation("regexp", "tests/data")
    i = a._conforming_index()
    assert i == -1

    genomepy.utils.rm_rf(gtf_file)
    genomepy.utils.rm_rf(bed_file)


def test__contig_conversion_dict():
    gtf_file = "tests/data/data.annotation.gtf"
    with open(gtf_file, "w") as f:
        f.write("\n")
    bed_file = "tests/data/data.annotation.bed"
    with open(bed_file, "w") as f:
        f.write("\n")
    genome_file = "tests/data/data.fa"
    with open(genome_file, "w") as f:
        f.write(">this matches chr2\n" ">that matches chr3\n")
    a = genomepy.annotation.Annotation("data", "tests")

    conversion_dict, duplicate_contigs = a._contig_conversion_dict(2)
    assert not duplicate_contigs
    assert conversion_dict["chr2"] == "this"
    assert conversion_dict["chr3"] == "that"

    # assume "matches" is the contig name
    ids, duplicate_contigs = a._contig_conversion_dict(1)
    assert "matches" in duplicate_contigs
    assert ids["matches"] == "that"

    genomepy.utils.rm_rf(gtf_file)
    genomepy.utils.rm_rf(bed_file)
    genomepy.utils.rm_rf(genome_file)


def test__conform_gtf():
    bed_file = "tests/data/data.annotation.bed"
    with open(bed_file, "w") as f:
        f.write("\n")

    # filter missing contigs
    conversion_dict = {"chr1": "chromosome_one"}
    gtf_file = "tests/data/data.annotation.gtf"
    with open(gtf_file, "w") as f:
        f.write("chr1\t1\t2\n")
        f.write("chr2\t1\t2\n")
    a = genomepy.annotation.Annotation("data", "tests")
    missing_contigs = a._conform_gtf(conversion_dict)
    assert missing_contigs == ["chr2"]
    with open(gtf_file) as f:
        lines = f.readlines()
    assert lines == ["chromosome_one\t1\t2\n"]

    # keep missing contigs
    with open(gtf_file, "w") as f:
        f.write("chr1\t1\t2\n")
        f.write("chr2\t1\t2\n")
    a = genomepy.annotation.Annotation("data", "tests")
    missing_contigs = a._conform_gtf(conversion_dict, filter_contigs=False)
    assert missing_contigs == ["chr2"]
    with open(gtf_file) as f:
        lines = f.readlines()
    assert lines == ["chromosome_one\t1\t2\n", "chr2\t1\t2\n"]

    genomepy.utils.rm_rf(gtf_file)
    genomepy.utils.rm_rf(bed_file)


def test_bed_from_gtf():
    bed_file = "tests/data/data.annotation.bed.gz"
    with genomepy.annotation._open(bed_file, "w") as f:
        f.write("this is the old bed\n")

    gtf_file = "tests/data/data.annotation.gtf.gz"
    with genomepy.annotation._open(gtf_file, "w") as f:
        f.write(
            "chr1\tgenomepy\texon\t1\t100\t.\t+\t.\t"
            'gene_id "ENSGP1234"; transcript_id "GP_1234.1";  gene_name "GP";\n'
        )
    a = genomepy.annotation.Annotation("data", "tests")

    # BED file before & after
    with genomepy.annotation._open(a.annotation_bed_file, "r") as f:
        lines = f.readlines()
    assert lines == ["this is the old bed\n"]
    a.bed_from_gtf()
    with genomepy.annotation._open(a.annotation_bed_file, "r") as f:
        lines = f.readlines()
    assert lines == ["chr1\t0\t100\tGP_1234.1\t0\t+\t100\t100\t0\t1\t100,\t0,\n"]
    assert len(lines) == 1

    # unzipping and rezipping correctly
    assert a.annotation_gtf_file.endswith(".gz")
    assert a.annotation_bed_file.endswith(".gz")

    genomepy.utils.rm_rf(gtf_file)
    genomepy.utils.rm_rf(bed_file)


def test_sanitize(caplog):
    tmp_dir = mkdtemp(dir="tests")
    genome = "testgenome"
    genomepy.utils.mkdir_p(os.path.join(tmp_dir, genome))

    bed_file = os.path.join(tmp_dir, genome, genome + ".annotation.bed")
    gtf_file = os.path.join(tmp_dir, genome, genome + ".annotation.gtf")
    genome_file = os.path.join(tmp_dir, genome, genome + ".fa")
    sizes_file = os.path.join(tmp_dir, genome, genome + ".fa.sizes")

    # no genome
    with open(bed_file, "w") as f:
        f.write("\n")
    with open(gtf_file, "w") as f:
        f.write("\n")
    a = genomepy.annotation.Annotation(genome, tmp_dir)
    a.sanitize()
    assert "A genome is required for sanitizing!" in caplog.text

    # conforming, filtering off
    with open(bed_file, "w") as f:
        f.write("chr1\t0\t100\n" "chr2\t0\t100\n")
    with open(gtf_file, "w") as f:
        f.write(
            "chr1\tgenomepy\texon\t1\t100\t.\t+\t.\t"
            'gene_id "ENSGP1234"; transcript_id "GP_1234.1";  gene_name "GP1";\n'
            "chr2\tgenomepy\texon\t1\t100\t.\t+\t.\t"
            'gene_id "ENSGP1235"; transcript_id "GP_1235.1";  gene_name "GP2";\n'
        )
    with open(genome_file, "w") as f:
        f.write(">chr1\n")
        f.write("ATCGATCG\n")
    genomepy.files.generate_fa_sizes(genome_file, sizes_file)
    a = genomepy.annotation.Annotation(genome, tmp_dir)
    a.sanitize(match_contigs=False, filter_contigs=False)

    # expect no changes
    with open(gtf_file) as f:
        lines = f.readlines()
    assert lines == [
        'chr1\tgenomepy\texon\t1\t100\t.\t+\t.\tgene_id "ENSGP1234"; transcript_id "GP_1234.1";  gene_name "GP1";\n',
        'chr2\tgenomepy\texon\t1\t100\t.\t+\t.\tgene_id "ENSGP1235"; transcript_id "GP_1235.1";  gene_name "GP2";\n',
    ]
    with open(bed_file) as f:
        lines = f.readlines()
    assert lines == ["chr1\t0\t100\n", "chr2\t0\t100\n"]
    metadata, _ = genomepy.files.read_readme(a.readme_file)
    assert metadata["sanitized annotation"] == "contigs match but not filtered"

    # conforming, filtering on
    a.sanitize(match_contigs=False, filter_contigs=True)

    # expect chr2 to have been removed
    with open(gtf_file) as f:
        lines = f.readlines()
    assert lines == [
        'chr1\tgenomepy\texon\t1\t100\t.\t+\t.\tgene_id "ENSGP1234"; transcript_id "GP_1234.1";  gene_name "GP1";\n',
    ]
    with open(bed_file) as f:
        lines = f.readlines()
    assert lines == [
        "chr1\t0\t100\tGP_1234.1\t0\t+\t100\t100\t0\t1\t100,\t0,\n",
    ]
    metadata, _ = genomepy.files.read_readme(a.readme_file)
    assert metadata["sanitized annotation"] == "contigs match and filtered"

    # not conforming, no fix possible
    genomepy.utils.rm_rf(os.path.join(tmp_dir, genome))
    genomepy.utils.mkdir_p(os.path.join(tmp_dir, genome))
    with open(bed_file, "w") as f:
        f.write("chr1\t0\t100\n")
    with open(gtf_file, "w") as f:
        f.write("\n")
    with open(genome_file, "w") as f:
        f.write(">not_chr1\n")
        f.write("ATCGATCG\n")
    genomepy.files.generate_fa_sizes(genome_file, sizes_file)
    a = genomepy.annotation.Annotation(genome, tmp_dir)
    a.sanitize()
    metadata, _ = genomepy.files.read_readme(a.readme_file)
    assert metadata["sanitized annotation"] == "not possible"

    # not conforming, filtering off
    genomepy.utils.rm_rf(os.path.join(tmp_dir, genome))
    genomepy.utils.mkdir_p(os.path.join(tmp_dir, genome))
    with open(bed_file, "w") as f:
        f.write("chr1\t0\t100\n" "chr2\t0\t100\n")
    with open(gtf_file, "w") as f:
        f.write(
            "chr1\tgenomepy\texon\t1\t100\t.\t+\t.\t"
            'gene_id "ENSGP1234"; transcript_id "GP_1234.1";  gene_name "GP1";\n'
            "chr2\tgenomepy\texon\t1\t100\t.\t+\t.\t"
            'gene_id "ENSGP1235"; transcript_id "GP_1235.1";  gene_name "GP2";\n'
        )
    with open(genome_file, "w") as f:
        f.write(">this matches chr1\n")
        f.write("ATCGATCG\n")
        f.write(">this2 matches chr1\n")
        f.write("ATCGATCG\n")
    genomepy.files.generate_fa_sizes(genome_file, sizes_file)
    a = genomepy.annotation.Annotation(genome, tmp_dir)
    a.sanitize(match_contigs=True, filter_contigs=False)
    assert "The genome contains duplicate contig names" in caplog.text
    with open(bed_file) as f:
        lines = f.readlines()
    assert lines == [
        "this2\t0\t100\tGP_1234.1\t0\t+\t100\t100\t0\t1\t100,\t0,\n",
        "chr2\t0\t100\tGP_1235.1\t0\t+\t100\t100\t0\t1\t100,\t0,\n",
    ]
    metadata, _ = genomepy.files.read_readme(a.readme_file)
    assert metadata["sanitized annotation"] == "contigs fixed but not filtered"

    # not conforming, filtering off
    genomepy.utils.rm_rf(os.path.join(tmp_dir, genome))
    genomepy.utils.mkdir_p(os.path.join(tmp_dir, genome))
    with open(bed_file, "w") as f:
        f.write("chr1\t0\t100\n" "chr2\t0\t100\n")
    with open(gtf_file, "w") as f:
        f.write(
            "chr1\tgenomepy\texon\t1\t100\t.\t+\t.\t"
            'gene_id "ENSGP1234"; transcript_id "GP_1234.1";  gene_name "GP1";\n'
            "chr2\tgenomepy\texon\t1\t100\t.\t+\t.\t"
            'gene_id "ENSGP1235"; transcript_id "GP_1235.1";  gene_name "GP2";\n'
        )
    with open(genome_file, "w") as f:
        f.write(">this matches chr1\n")
        f.write("ATCGATCG\n")
        f.write(">this2 matches chr1\n")
        f.write("ATCGATCG\n")
    genomepy.files.generate_fa_sizes(genome_file, sizes_file)
    a = genomepy.annotation.Annotation(genome, tmp_dir)
    a.sanitize(match_contigs=True, filter_contigs=True)
    assert "The genome contains duplicate contig names" in caplog.text
    with open(bed_file) as f:
        lines = f.readlines()
    assert lines == ["this2\t0\t100\tGP_1234.1\t0\t+\t100\t100\t0\t1\t100,\t0,\n"]
    metadata, _ = genomepy.files.read_readme(a.readme_file)
    assert metadata["sanitized annotation"] == "contigs fixed and filtered"

    genomepy.utils.rm_rf(tmp_dir)


def test_query_mygene():
    result = query_mygene(["ENST00000449992"], 9606, "symbol")
    assert result.loc["ENST00000449992", "symbol"] == "TP63"
