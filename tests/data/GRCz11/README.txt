name: GRCz11
provider: Ensembl
original name: GRCz11
original filename: na
assembly_accession: GCA_000002035.4
tax_id: 7955
mask: na
genome url: na
annotation url: ftp://ftp.ensembl.org/pub/release-104/gtf/danio_rerio/Danio_rerio.GRCz11.104.gtf.gz
sanitized annotation: na
genomepy version: na
date: na

annotation files were subset with this code
```
import genomepy
import csv

a = genomepy.Annotation("GRCz11", "tests/data")
bed = a.bed
bed = bed[~bed.chrom.str.contains("K")]  # drop scaffolds
bed = bed.groupby(["chrom"]).head(5)  # subset to 5 genes per chromosome
genomepy.annotations.utils.write_annot(
    bed,
    "tests/data/GRCz11/GRCz11.annotation.bed",
)

gtf = a.gtf
gtf = gtf[~gtf.seqname.str.contains("K")]  # drop scaffolds
gtf = gtf.groupby(["seqname"]).head(5)  # subset to 5 genes per chromosome
genomepy.annotations.utils.write_annot(
    gtf,
    "tests/data/GRCz11/GRCz11.annotation.gtf",
)
```
