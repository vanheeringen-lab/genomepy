# genomepy

[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io)
[![PyPI version](https://badge.fury.io/py/genomepy.svg)](https://badge.fury.io/py/genomepy)
[![Build Status](https://travis-ci.org/simonvh/genomepy.svg?branch=master)](https://travis-ci.org/simonvh/genomepy)
[![Code Health](https://landscape.io/github/simonvh/genomepy/master/landscape.svg?style=flat)](https://landscape.io/github/simonvh/genomepy/master)

 [![status](http://joss.theoj.org/papers/df434a15edd00c8c2f4076668575d1cd/status.svg)](http://joss.theoj.org/papers/df434a15edd00c8c2f4076668575d1cd)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.831969.svg)](https://doi.org/10.5281/zenodo.831969)

Easily install and use genomes in Python and elsewhere!

The goal is to have a _simple_ and _straightforward_ way to download and use genomic sequences. 
Currently, genomepy supports UCSC, Ensembl and NCBI. 

## Installation

Genomepy works with Python 2.7 and Python 3.4+. 
You can install it via [bioconda](https://bioconda.github.io/):

```
$ conda install genomepy
``` 

Or via pip:

```
$ pip install genomepy
```

If you install via pip, you will have to install some dependencies, only if 
you want to use the annotation download feature. 
You will have to install the following
utilities and make sure they are in your PATH:

* `genePredToBed`
* `genePredToGtf`
* `bedToGenePred`
* `gtfToGenePred`

You can find the binaries [here](http://hgdownload.cse.ucsc.edu/admin/exe/).

## Configuration

By default genomes will be saved in `~/.local/share/genomes`. 
This default can be changed by creating a configuration file called `~/.config/genomepy/genomepy.yaml`. 
For instance, to set the default genome directory to `/data/genomes`, edit `~/.config/genomepy/genomepy.yaml` and add the following line:

```
genome_dir: /data/genomes
```

The genome directory can also be explicitly specified in both the Python API as well as on the command-line.

## Usage

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


#### List available providers

```
$ genomepy providers
Ensembl
UCSC
NCBI
```

#### List available genomes

You can constrain the genome list by using the `-p` option to search only a 
specific provider. 

```
$ genomepy genomes -p UCSC
UCSC	hg38	Human Dec. 2013 (GRCh38/hg38) Genome at UCSC
UCSC	hg19	Human Feb. 2009 (GRCh37/hg19) Genome at UCSC
UCSC	hg18	Human Mar. 2006 (NCBI36/hg18) Genome at UCSC
...
UCSC	danRer4	Zebrafish Mar. 2006 (Zv6/danRer4) Genome at UCSC
UCSC	danRer3	Zebrafish May 2005 (Zv5/danRer3) Genome at UCSC
```

#### Search for a genome.

```
$ genomepy search Xenopus
NCBI	Xenopus_tropicalis_v9.1	Xenopus tropicalis; DOE Joint Genome Institute
NCBI	ViralProj30173	Xenopus laevis endogenous retrovirus Xen1; 
NCBI	Xenopus_laevis_v2	Xenopus laevis; International Xenopus Sequencing Consortium
NCBI	v4.2	Xenopus tropicalis; DOE Joint Genome Institute
NCBI	Xtropicalis_v7	Xenopus tropicalis; DOE Joint Genome Institute
Ensembl	JGI 4.2	Xenopus
```

Only search a specific provider:

```
$ genomepy search tropicalis -p UCSC
UCSC	xenTro7	X. tropicalis Sep. 2012 (JGI 7.0/xenTro7) Genome at UCSC
UCSC	xenTro3	X. tropicalis Nov. 2009 (JGI 4.2/xenTro3) Genome at UCSC
UCSC	xenTro2	X. tropicalis Aug. 2005 (JGI 4.1/xenTro2) Genome at UCSC
UCSC	xenTro1	X. tropicalis Oct. 2004 (JGI 3.0/xenTro1) Genome at UCSC
```

Note that searching doesn't work flawlessly, so try a few variations if 
you don't get any results. 
Search is case-insensitive.


#### Install a genome.

The most important command. The most simple form:

```
$ genomepy  install hg38 UCSC 
downloading...
done...
name: hg38
fasta: /data/genomes/hg38/hg38.fa
```

Here, genomes are downloaded to the directory specified in the config file. 
To choose a different directory, use the `-g` option.

```
$ genomepy install sacCer3 UCSC -g ~/genomes/
downloading from http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz...
done...
name: sacCer3
local name: sacCer3
fasta: /home/simon/genomes/sacCer3/sacCer3.fa
```

You can use a regular expression to filter for matching sequences 
(or non-matching sequences by using the `--no-match` option). For instance, 
the following command downloads hg38 and saves only the major chromosomes:

```
$ genomepy  install hg38 UCSC -r 'chr[0-9XY]+$'
downloading from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz...
done...
name: hg38
local name: hg38
fasta: /data/genomes/hg38/hg38.fa
$ grep ">" /data/genomes/hg38/hg38.fa
>chr1
>chr10
>chr11
>chr12
>chr13
>chr14
>chr15
>chr16
>chr17
>chr18
>chr19
>chr2
>chr20
>chr21
>chr22
>chr3
>chr4
>chr5
>chr6
>chr7
>chr8
>chr9
>chrX
>chrY
```

By default, sequences are soft-masked. Use `-m hard` for hard masking.

The chromosome sizes are saved in file called `<genome_name>.fa.sizes`.

For genomes from UCSC and Ensembl, you can choose to download gene annotation
files with the `--annotation` option. 
These will be saved in BED and GTF format. 

```
$ genomepy  install hg38 UCSC --annotation
```

Finally, in the spirit of reproducibility all selected options are stored in a `README.txt`. 
This includes the original name and download location. 

#### Local cache. 

Note that the first time you run `genomepy search` or `list` the command will take a long time
as the genome lists have to be downloaded. 
The lists are cached locally, which will save time later. The cached files are stored in 
`~/.cache/genomepy` and expire after 7 days. You can also delete this directory to clean the 
cache.

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

The `genomepy.genome()` method returns a Genome object. This has all the
functionality of a `pyfaidx.Fasta` object, 
see the [documentation](https://github.com/mdshw5/pyfaidx) for more examples on how to use this.

## Known issues

There might be issues with specific genome sequences.
Sadly, not everything (naming, structure, filenames) is always consistent on the provider end. 
Let me know if you encounter issues with certain downloads.

## Todo

* Linking genomes to NCBI taxonomy ID
* Optionally: automatic indexing for aligners (such as bwa)
* Optionally: Ensembl bacteria (although there might be better options specifically for bacterial sequences)

## Citation

If you use genomepy in your research, please consider citing it: [10.21105/joss.00320](http://dx.doi.org/10.21105/joss.00320).


## Getting help

If you want to report a bug or issue, or have problems with installing or running the software please create [a new issue](https://github.com/simonvh/genomepy/issues). This is the preferred way of getting support. Alternatively, you can [mail me](mailto:simon.vanheeringen@gmail.com).

## Contributing

Contributions welcome! Send me a pull request or get in [touch](mailto:simon.vanheeringen@gmail.com).

## License

This module is licensed under the terms of the [MIT license](https://opensource.org/licenses/MIT).
