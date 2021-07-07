genome: GRCz11
provider: NCBI

GTF creation:
tail -n 4 GRCz11.annotation.gtf > tests/data/sanitize/sanitize.annotation.gtf
head -n 4 GRCz11.annotation.gtf >> tests/data/sanitize/sanitize.annotation.gtf

BED creation:
genomepy.annotation.utils.generate_annot("tests/data/sanitize/sanitize.annotation.gtf", "tests/data/sanitize/sanitize.annotation.bed")

FASTA creation:
cat GRCz11.fa | grep "" > tests/data/sanitize/sanitize.fa
echo ACTG >> tests/data/sanitize/sanitize.fa
