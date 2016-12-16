# genomepy

[![PyPI version](https://badge.fury.io/py/genomepy.svg)](https://badge.fury.io/py/genomepy)

Easily install and use genomes in Python and elsewhere!

## Installation

Via pip, for now.

```
$ pip install genomepy
```

## Usage

```python
>>> import genomepy
>>> for row in genomepy.search("human"):
...     print "\t".join(row)
...
UCSC	hg38	Human Dec. 2013 (GRCh38/hg38) Genome at UCSC
UCSC	hg19	Human Feb. 2009 (GRCh37/hg19) Genome at UCSC
UCSC	hg18	Human Mar. 2006 (NCBI36/hg18) Genome at UCSC
UCSC	hg17	Human May 2004 (NCBI35/hg17) Genome at UCSC
UCSC	hg16	Human July 2003 (NCBI34/hg16) Genome at UCSC
Ensembl	bacteria_102_collection_core_34_87_1	Brucella melitensis (GCA_000988815)
Ensembl	bacteria_94_collection_core_34_87_1	Brucella suis (GCA_000875695)
Ensembl	bacteria_131_collection_core_34_87_1	Candidatus Paraburkholderia schumannianae
Ensembl	homo_sapiens_core_86_38	Human
Ensembl	pediculus_humanus_core_34_87_2	Pediculus humanus
>>> genomepy.install_genome("hg38", "UCSC", "/data/genomes")
downloading...
done...
name: hg38
fasta: /data/genomes/hg38/hg38.fa
>>> g = genomepy.genome("hg38", "/data/genomes")
>>> g["chr6"][166502000:166503000]
tgtatggtccctagaggggccagagtcacagagatggaaagtggatggcgggtgccgggggctggggagctactgtgcagggggacagagctttagttctgcaagatgaaacagttctggagatggacggtggggatgggggcccagcaatgggaacgtgcttaatgccactgaactgggcacttaaacgtggtgaaaactgtaaaagtcatgtgtatttttctacaattaaaaaaaATCTGCCACAGAGTTAAAAAAATAACCACTATTTTCTGGAAATGGGAAGGAAAAGTTACAGCATGTAATTAAGATGACAATTTATAATGAACAAGGCAAATCTTTTCATCTTTGCCTTTTGGGCATATTCAATCTTTGCCCAGAATTAAGCACCTTTCAAGATTAATTCTCTAATAATTCTAGTTGAACAACACAACCTTTTCCTTCAAGCTTGCAATTAAATAAGGCTATTTTTAGCTGTAAGGATCACGCTGACCTTCAGGAGCAATGAGAACCGGCACTCCCGGCCTGAGTGGATGCACGGGGAGTGTGTCTAACACACAGGCGTCAACAGCCAGGGCCGCACGAGGAGGAGGAGTGGCAACGTCCACACAGACTCACAACACGGCACTCCGACTTGGAGGGTAATTAATACCAGGTTAACTTCTGGGATGACCTTGGCAACGACCCAAGGTGACAGGCCAGGCTCTGCAATCACCTCCCAATTAAGGAGAGGCGAAAGGGGACTCCCAGGGCTCAGAGCACCACGGGGTTCTAGGTCAGACCCACTTTGAAATGGAAATCTGGCCTTGTGCTGCTGCTCTTGTGGGGAGACAGCAGCTGCGGAGGCTGCTCTCTTCATGGGATTACTCTGGATAAAGTCTTTTTTGATTCTACgttgagcatcccttatctgaaatgcctgaaaccggaagtgtttaggatttggggattttgcaatatttacttatatataatgagatatcttggagatgggccacaa
```

The `genomepy.genome()` method returns a `pyfaidx.Fasta` object, 
see the [documentation](https://github.com/mdshw5/pyfaidx) for more examples on how to use this.

### Command line 

```
$ genomepy

Usage: genomepy [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  genomes    list available genomes
  install    install genome
  providers  list available providers
  search     search for genomes
```

List available genomes.

```
$ genomepy genomes -p UCSC
UCSC	hg38	Human Dec. 2013 (GRCh38/hg38) Genome at UCSC
UCSC	hg19	Human Feb. 2009 (GRCh37/hg19) Genome at UCSC
UCSC	hg18	Human Mar. 2006 (NCBI36/hg18) Genome at UCSC
...
UCSC	danRer4	Zebrafish Mar. 2006 (Zv6/danRer4) Genome at UCSC
UCSC	danRer3	Zebrafish May 2005 (Zv5/danRer3) Genome at UCSC
```

Install a genome.

```
$ genomepy  install hg38 UCSC /data/genomes/
downloading...
done...
name: hg38
fasta: /data/genomes/hg38/hg38.fa
```

## Todo

* Tests!
* Ensembl bacteria
* Automatic indexing (such as bwa)
* Caching of UCSC/Ensembl genome listings
* Configurable default genome installation directory

## License

This module is licensed under the terms of the [MIT license](https://opensource.org/licenses/MIT).
