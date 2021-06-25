name: sacCer3
provider: UCSC
original name: sacCer3
original filename: chromFa.tar.gz
assembly_accession: GCA_000146045.2
tax_id: 559292
mask: soft
genome url: http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz
annotation url: http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/genes/sacCer3.ensGene.gtf.gz
sanitized annotation: contigs match
genomepy version: 0.9.3
date: na
regex: 'alt' (inverted match)

No contigs were removed.

files were subset with this code
```
import genomepy
import csv

a = genomepy.Annotation("sacCer3", "tests/data")
bed = a.bed
bed = bed[bed.chrom == "chrIV"]
genes = list(set(bed.name))[0:50]
bed = bed[bed.name.isin(genes)]
bed.to_csv(
    "tests/data/sacCer3/sacCer3.annotation.bed",
    sep="\t",
    header=None,
    index=None,
    quoting=csv.QUOTE_NONE,
)

gtf = a.gtf
gtf = gtf[gtf.seqname == "chrIV"]
gtf2 = gtf[gtf.attribute.str.contains('|'.join(genes))]
gtf.to_csv(
    "tests/data/sacCer3/sacCer3.annotation.gtf",
    sep="\t",
    header=None,
    index=None,
    quoting=csv.QUOTE_NONE,
)

genomepy.files.filter_fasta("tests/data/sacCer3/sacCer3.fa", "chrIV", outfa="tests/data/sacCer3/sacCer3.fa")
```
