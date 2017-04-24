# genomepy

[![PyPI version](https://badge.fury.io/py/genomepy.svg)](https://badge.fury.io/py/genomepy)
[![Build Status](https://travis-ci.org/simonvh/genomepy.svg?branch=master)](https://travis-ci.org/simonvh/genomepy)
[![Code Health](https://landscape.io/github/simonvh/genomepy/master/landscape.svg?style=flat)](https://landscape.io/github/simonvh/genomepy/master)

Easily install and use genomes in Python and elsewhere!

The goal is to have a _simple_ and _straightforward_ way to download and use genomic sequences. 
Currently, genomepy supports UCSC, Ensembl and NCBI. 

## Installation

Via pip, for now.

```
$ pip install genomepy
```

## Configuration

By default genomes will be saved in `~/.local/share/genomes`. 
This default can be changed by creating a configuration file called `~/.config/genomepy/genomepy.yaml`. 
For instance, to set the default genome directory to `/data/genomes`, edit `~/.config/genomepy/genomepy.yaml` and add the following line:

```
genome_dir: /data/genomes
```

The genome directory can also be explicitly specified in both the Python API as well as on the command-line.

## Usage

### From Python

```python
>>> import genomepy
>>> for row in genomepy.search("GRCh38"):
...     print("\t".join(row))
...
UCSC	hg38	Human Dec. 2013 (GRCh38/hg38) Genome at UCSC
NCBI	GRCh38.p10	Homo sapiens; Genome Reference Consortium
NCBI	GRCh38	Homo sapiens; Genome Reference Consortium
NCBI	GRCh38.p1	Homo sapiens; Genome Reference Consortium
NCBI	GRCh38.p2	Homo sapiens; Genome Reference Consortium
NCBI	GRCh38.p3	Homo sapiens; Genome Reference Consortium
NCBI	GRCh38.p4	Homo sapiens; Genome Reference Consortium
NCBI	GRCh38.p5	Homo sapiens; Genome Reference Consortium
NCBI	GRCh38.p6	Homo sapiens; Genome Reference Consortium
NCBI	GRCh38.p7	Homo sapiens; Genome Reference Consortium
NCBI	GRCh38.p8	Homo sapiens; Genome Reference Consortium
NCBI	GRCh38.p9	Homo sapiens; Genome Reference Consortium
Ensembl	GRCh38.p10	Human
>>> genomepy.install_genome("hg38", "UCSC", genome_dir="/data/genomes")
downloading...
done...
name: hg38
fasta: /data/genomes/hg38/hg38.fa
>>> g = genomepy.genome("hg38", genome_dir="/data/genomes")
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
$ genomepy  install hg38 UCSC -g /data/genomes
downloading...
done...
name: hg38
fasta: /data/genomes/hg38/hg38.fa
```

## Known issues

There might be issues with specific genome sequences.
Sadly, not everything (naming, structure, filenames) is always consistent on the provider end. 
Let me know if you encounter issues with certain downloads.

## Todo

* More tests!
* Caching of genome listings
* Automatic indexing (such as bwa)
* Ensembl bacteria

## Contributing

Contributions welcome! Send me a pull request or get in [touch](mailto:simon.vanheeringen@gmail.com).

## License

This module is licensed under the terms of the [MIT license](https://opensource.org/licenses/MIT).
