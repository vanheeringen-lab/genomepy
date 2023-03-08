# genomepy: genes and genomes at your fingertips

[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/genomepy/badges/downloads.svg)](https://anaconda.org/bioconda/genomepy)
[![PyPI version](https://badge.fury.io/py/genomepy.svg)](https://badge.fury.io/py/genomepy)
[![GitHub stars](https://badgen.net/github/stars/vanheeringen-lab/genomepy)](https://GitHub.com/vanheeringen-lab/genomepy/stargazers/)

[![Build Status](https://app.travis-ci.com/vanheeringen-lab/genomepy.svg?branch=master)](https://app.travis-ci.com/github/vanheeringen-lab/genomepy/branches)
[![Maintainability](https://api.codeclimate.com/v1/badges/c4476820f1d21a3e0569/maintainability)](https://codeclimate.com/github/vanheeringen-lab/genomepy/maintainability)
[![Test Coverage](https://api.codeclimate.com/v1/badges/c4476820f1d21a3e0569/test_coverage)](https://codeclimate.com/github/vanheeringen-lab/genomepy/test_coverage)

[![bioinformatics](https://img.shields.io/badge/DOI-10.1093%2Fbioinformatics%2Fbtad119-%23167da4)](https://doi.org/10.1093/bioinformatics/btad119)
[![zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.1010458.svg)](https://doi.org/10.5281/zenodo.1010458)
<!--- [![arXiv](https://img.shields.io/badge/arXiv-10.48550/arXiv.2209.00842-b31b1b.svg)](https://doi.org/10.48550/arXiv.2209.00842) [![status](http://joss.theoj.org/papers/df434a15edd00c8c2f4076668575d1cd/status.svg)](http://joss.theoj.org/papers/df434a15edd00c8c2f4076668575d1cd) -->

genomepy is designed to provide a _simple_ and _straightforward_ way to download and use genomic data.
This includes (1) searching available data, 
(2) showing the available metadata,
(3) automatically downloading, preprocessing and matching data and
(4) generating optional aligner indexes.
All with sensible, yet controllable defaults.
Currently, genomepy supports Ensembl, UCSC, NCBI and GENCODE. 

[![asciicast](https://asciinema.org/a/eZttBuf5ly0AnjFVBiEIybbjS.png)](https://asciinema.org/a/eZttBuf5ly0AnjFVBiEIybbjS)

**Pssst, hey there!** Is genomepy not doing what you want? Does it fail? Is it clunky?
Is the documentation unclear? Have any other ideas on how to improve it?
Don't be shy and [let us know](https://github.com/vanheeringen-lab/genomepy/issues)!

## Table of Contents
1.  [Installation](#installation)
2.  [Quick usage](#quick-usage)
3.  [Command line](#command-line-interface)
    1. [Search](#search-genomes--gene-annotations)
    2. [Inspect](#inspect-gene-annotations)
    3. [Install](#install-a-genome--gene-annotation)
    4. [Plugins](#plugins)
4.  [Python](#python-api)
5.  [Frequently Asked Questions](#frequently-asked-questions)
6.  [Getting help](#getting-help)
7.  [Contributing](#contributing)
8.  [Citation](#citation)
9.  [License](#license)

## Installation

genomepy requires Python 3.7+

You can install genomepy via [bioconda](https://bioconda.github.io/), pip or git.

#### Bioconda

```bash
$ conda install -c bioconda genomepy
``` 

#### Pip

```bash
$ pip install genomepy
```

With the Pip installation, you will have to install additional dependencies, and make them available in your PATH.

To read/write bgzipped genomes you will have to install `pysam`.

If you want to use gene annotation features, you will have to install the following utilities:

* `genePredToBed`
* `genePredToGtf`
* `bedToGenePred`
* `gtfToGenePred`
* `gff3ToGenePred`

You can find the binaries [here](http://hgdownload.cse.ucsc.edu/admin/exe/).

#### Git

```bash
$ git clone https://github.com/vanheeringen-lab/genomepy.git
$ conda env create -n genomepy -f genomepy/environment.yml
$ conda activate genomepy
$ pip install -e genomepy
```

## Quick usage

1.  Find your genome: `$ genomepy search zebrafish`

  Console output:
  ```bash
  name      provider    accession           tax_id     annotation    species        other_info
  GRCz11    Ensembl     GCA_000002035.4     7955       ✓             Danio rerio    2017-08-Ensembl/2018-04
   ^
   Use name for genomepy install
  ```

2.  Install your genome (with gene annotation): `$ genomepy install --annotation GRCz11 --provider ensembl `

The default genomes directory: `~/.local/share/genomes/`

## Command line interface

All commands come with a short explanation when appended with `-h`/`--help`.

```bash
$ genomepy --help
Usage: genomepy [OPTIONS] COMMAND [ARGS]...

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  annotation  show 1st lines of each annotation
  clean       remove provider data
  config      manage configuration
  genomes     list available genomes
  install     install a genome & run active plugins
  plugin      manage plugins
  providers   list available providers
  search      search for genomes
```

### Search genomes & gene annotations

Let's say we want to download a *Xenopus tropicalis* genome & gene annotation. 
First, lets find out what's out there!

You can search by name, taxonomy ID or assembly accession ID.
Additionally, you can limit the search result to one provider with `-p`/`--provider`.
Furthermore, you can get the absolute `--size` of each genome (this option slows down the search).

```bash
$ genomepy search xenopus tro
name                       provider    accession           tax_id    annotation     species               other_info
                                                                      n r e k
Xenopus_tropicalis_v9.1    Ensembl     GCA_000004195.3       8364        ✓          Xenopus tropicalis    2019-04-Ensembl/2019-12
xenTro1                    UCSC        na                    8364     ✗ ✗ ✗ ✗       Xenopus tropicalis    Oct. 2004 (JGI 3.0/xenTro1)
xenTro2                    UCSC        na                    8364     ✗ ✓ ✓ ✗       Xenopus tropicalis    Aug. 2005 (JGI 4.1/xenTro2)
xenTro3                    UCSC        GCA_000004195.1       8364     ✗ ✓ ✓ ✗       Xenopus tropicalis    Nov. 2009 (JGI 4.2/xenTro3)
xenTro7                    UCSC        GCA_000004195.2       8364     ✓ ✓ ✗ ✗       Xenopus tropicalis    Sep. 2012 (JGI 7.0/xenTro7)
xenTro9                    UCSC        GCA_000004195.3       8364     ✓ ✓ ✓ ✗       Xenopus tropicalis    Jul. 2016 (Xenopus_tropicalis_v9.1/xenTro9)
Xtropicalis_v7             NCBI        GCF_000004195.2       8364        ✓          Xenopus tropicalis    DOE Joint Genome Institute
Xenopus_tropicalis_v9.1    NCBI        GCF_000004195.3       8364        ✓          Xenopus tropicalis    DOE Joint Genome Institute
UCB_Xtro_10.0              NCBI        GCF_000004195.4       8364        ✓          Xenopus tropicalis    University of California, Berkeley
ASM1336827v1               NCBI        GCA_013368275.1       8364        ✗          Xenopus tropicalis    Southern University of Science and Technology
 ^
 Use name for genomepy install
```

### Inspect gene annotations

Let's say we want to download the *Xenopus tropicalis* genome & gene annotation from UCSC.

Since we are interested in the gene annotation as well, we should check which gene annotation suits our needs.
As you can see in the search results, UCSC has several gene annotations for us to choose from.
In the search results, `n r e k ` denotes which UCSC annotations are available. 
These stand for **n**cbiRefSeq, **r**efGene, **e**nsGene and **k**nownGene, respectively.

We can quickly inspect these with the `genomepy annotation` command:

```bash
$ genomepy annotation xenTro9 -p ucsc
12:04:41 | INFO | UCSC ncbiRefSeq
chr1    genomepy        transcript      133270  152620  .       -       .       gene_id "LOC100490505"; transcript_id "XM_012956089.1";  gene_name "LOC100490505";
chr1    genomepy        exon    133270  134186  .       -       .       gene_id "LOC100490505"; transcript_id "XM_012956089.1"; exon_number "1"; exon_id "XM_012956089.1.1"; gene_name "LOC100490505";
12:04:45 | INFO | UCSC refGene
chr1    genomepy        transcript      193109390       193134311       .       +       .       gene_id "pias2"; transcript_id "NM_001078987";  gene_name "pias2";
chr1    genomepy        exon    193109390       193109458       .       +       .       gene_id "pias2"; transcript_id "NM_001078987"; exon_number "1"; exon_id "NM_001078987.1"; gene_name "pias2";
12:04:49 | INFO | UCSC ensGene
chr1    genomepy        transcript      133270  152620  .       -       .       gene_id "ENSXETG00000030302.2"; transcript_id "ENSXETT00000061673.2";  gene_name "ENSXETG00000030302.2";
chr1    genomepy        exon    133270  134186  .       -       .       gene_id "ENSXETG00000030302.2"; transcript_id "ENSXETT00000061673.2"; exon_number "1"; exon_id "ENSXETT00000061673.2.1"; gene_name "ENSXETG00000030302.2";
```

Here we can see that the `refGene` annotation has actual HGNC gene names, so lets go with this annotation.
This differs between assemblies, so be sure to check!

### Install a genome & gene annotation

Copy the name returned by the search function to install.

```bash
$ genomepy install xenTro9
```

You can choose to download gene annotation files with the `-a`/`--annotation` option.

```bash
$ genomepy install xenTro9 --annotation
```

For UCSC we can also select the annotation type.
See `genomepy install --help` for all provider specific options.

```bash
$ genomepy install xenTro9 --UCSC-annotation refGene
```

Since we did not specify the provider here, genomepy will use the first provider with `xenTro9`.
You can specify a provider by name with `-p`/`--provider`:

```bash
$ genomepy install xenTro9 -p UCSC
Downloading genome from http://hgdownload.soe.ucsc.edu/goldenPath/xenTro9/bigZips/xenTro9.fa.gz...
Genome download successful, starting post processing...

name: xenTro9
local name: xenTro9
fasta: ~/.local/share/genomes/xenTro9/xenTro9.fa
```

Next, the genome is downloaded to the directory specified in the config file (by default `~/.local/share/genomes`).
To choose a different directory, use the `-g`/`--genomes_dir` option:

```bash
$ genomepy install sacCer3 -p UCSC -g /path/to/my/genomes
Downloading genome from http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz...
Genome download successful, starting post processing...

name: sacCer3
local name: sacCer3
fasta: /path/to/my/genomes/sacCer3/sacCer3.fa
```

#### Regex, masking & compression

You can use a regular expression to filter for matching sequences
(or non-matching sequences by using the `-n`/`--no-match` option).
For instance, the following command downloads hg38 and saves only the major chromosomes:

```bash
$ genomepy install hg38 -p UCSC -r 'chr[0-9XY]+$'
Downloading genome from from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz...
Genome download successful, starting post processing...

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

By default, genome sequences are soft-masked (ACgtN). 
Use `-m hard` for hard masking (ACNNN), or `-m none` for no masking (ACGTN).

```bash
$ genomepy install hg38 --mask hard
```

If you wish to conserve space, you can tell genomepy to compress the downloaded data by passing the `-b`/`--bgzip` option.
See [Configuration](#compression) for details.

```bash
$ genomepy install hg38 --bgzip
```

#### Other providers (any URL/local files)

To use assemblies not on NCBI, UCSC, Ensembl or GENCODE, you can give a URL instead of a name, together with `--provider URL`.
Similarly, if you have a local FASTA file, you can install this using the filepath, together with `--provider Local`:

```bash
$ genomepy install -p url https://research.nhgri.nih.gov/hydra/download/assembly/\Hm105_Dovetail_Assembly_1.0.fa.gz
```

This will install the genome under the filename of the URL/filepath, but can be changed with the `-l`/`--localname` option.

If you add the `--annotation` flag, genomepy will search the (remote) directory for an annotation file as well.
Should this fail, you can also add a URL to the annotation with `--URL-to-annotation` with the `URL` provider, 
or a filepath with `--Local-path-to-annotation` with the `Local` provider:

```bash
$ genomepy install -p local /path/to/genome.fa --Local-path-to-annotation /path/to/gene.annotation.gtf
```

#### Reproducibility

All selected options are stored in a `README.txt`.
This includes the original name, download location and other genomepy operations (such as regex filtering and time).

### Plugins

Plugins are optional steps that are executed after installing an assembly with `genomepy install`.
If you already installed an assembly, you can activate a plugin and rerun the install command.
This will not overwrite your local files, unless you use the `--force` option.

Check which plugins are enabled with `genomepy plugin list`.

#### Genome blacklists

For some model organisms, genomepy can download a genome blacklist (generated by the [Kundaje lab](https://www.nature.com/articles/s41598-019-45839-z)).
Blacklists are only available for these model organisms when downloaded from UCSC, and for the human and mouse genomes.

Enable the blacklist plugin to use it:

```bash
$ genomepy plugin enable blacklist
Enabled plugins: blacklist
```

#### Aligner indexes

You can also create aligner indexes for several widely used aligners.
Currently, genomepy supports:

* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [BWA](http://bio-bwa.sourceforge.net/)
* [GMAP](http://research-pub.gene.com/gmap/)
* [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)
* [Minimap2](https://github.com/lh3/minimap2)
* [STAR](https://github.com/alexdobin/STAR)

These programs are not installed by genomepy and need to be installed separately for the indexing to work.
The easiest way to do so is with conda, e.g.: `conda install -c bioconda bwa star`

Splice-aware indexing (required for e.g. RNA-seq) can be performed by STAR and Hisat2.
This will be done automatically if the gene annotation was downloaded as well.
Finally, STAR can further improve mapping to (novel) splice junctions by indexing again (see 2-pass mapping mode in the STAR manual).
The second pass is not supported by genomepy.

You can configure the index creation with `genomepy plugin enable`, e.g.:

```bash
$ genomepy plugin enable bwa star
Enabled plugins: blacklist, bwa, star
```

You can pass the number of threads to use for aligner index creation with `genomepy install --threads` (default is 8).

### Configuration

All defaults can be overwritten on the command line and in Python.
However, you can create & edit the config file to change the default settings:

```bash
$ genomepy config generate
Created config file ~/.config/genomepy/genomepy.yaml
```

#### Genome location

By default, genomes will be saved in `~/.local/share/genomes`.

To set the default genome directory, to `/data/genomes` for instance,
edit `~/.config/genomepy/genomepy.yaml` and change the following line:

```yaml
genomes_dir: /data/genomes
```

#### Compression

Genome FASTA files can be stored using bgzip compression.
This means that the FASTA files will take up less space on disk.
Set the following line to your config file:

```yaml
bgzip: True
```

Most tools are able to use bgzip-compressed genome files.
One notable exception is `bedtools getfasta`.
As an alternative, you can use the `faidx` command-line script from [pyfaidx](https://github.com/mdshw5/pyfaidx)
which comes installed with genomepy.

### List available providers

```bash
$ genomepy providers
GENCODE
Ensembl
UCSC
NCBI
Local
URL
```

### List available genomes

You can constrain the genome list by using the `-p`/`--provider` option to search only a specific provider.
Additionally, you can get the absolute `--size` of each genome (this option slows down the search).

```bash
$ genomepy genomes -p UCSC
name                    provider    accession          tax_id     annotation     species                                     other_info
                                                                   n r e k
ailMel1                 UCSC        GCF_000004335.2      9646      ✓ ✗ ✓ ✗       Ailuropoda melanoleuca                      Dec. 2009 (BGI-Shenzhen 1.0/ailMel1)
allMis1                 UCSC        GCA_000281125.1      8496      ✗ ✓ ✗ ✗       Alligator mississippiensis                  Aug. 2012 (allMis0.2/allMis1)
anoCar1                 UCSC        na                  28377      ✗ ✗ ✓ ✗       Anolis carolinensis                         Feb. 2007 (Broad/anoCar1)
```

### Local cache.

Note that the first time you run `genomepy search` or `list` the command will take a while as the genome lists have to be downloaded.
The lists are cached locally, which will save time later.
The cached files are stored in `~/.cache/genomepy` and expire after 7 days (so they stay up to date).
This expiration time can be changed in the config file.
You can also delete this directory to clean the cache using `genomepy clean`.

## Python API

Check out our [Python API documentation here](https://vanheeringen-lab.github.io/genomepy/content/api_core.html)

```
>>> import genomepy
>>> for row in genomepy.search("GRCh38"):
...    print(row)
...    
['GRCh38.p13', 'Ensembl', 'GCA_000001405.28', 9606, True, 'Homo sapiens', '2014-01-Ensembl/2021-03']
['hg38', 'UCSC', 'GCA_000001405.15', 9606, [True, True, False, True], 'Homo sapiens', 'Dec. 2013 (GRCh38/hg38)']
['GRCh38', 'NCBI', 'GCF_000001405.26', 9606, True, 'Homo sapiens', 'Genome Reference Consortium']
['GRCh38.p1', 'NCBI', 'GCF_000001405.27', 9606, True, 'Homo sapiens', 'Genome Reference Consortium']
['GRCh38.p2', 'NCBI', 'GCF_000001405.28', 9606, True, 'Homo sapiens', 'Genome Reference Consortium']
['GRCh38.p3', 'NCBI', 'GCF_000001405.29', 9606, True, 'Homo sapiens', 'Genome Reference Consortium']

>>> genomepy.install_genome("hg38", annotation=True, provider="UCSC", genomes_dir="./data/genomes")
Downloading genome from UCSC. Target URL: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz...
Genome download successful, starting post processing...
name: hg38
local name: hg38
fasta: ./data/genomes/hg38/hg38.fa
Downloading the ncbiRefSeq annotation from the UCSC MySQL database.
Annotation download successful

>>> a = genomepy.Annotation("hg38", genomes_dir="./data/genomes")
>>> a.named_gtf.head(3)
          seqname  ...                                          attribute
gene_name          ...                                                   
DDX11L1      chr1  ...  gene_id "DDX11L1"; transcript_id "NR_046018.2"...
DDX11L1      chr1  ...  gene_id "DDX11L1"; transcript_id "NR_046018.2"...
DDX11L1      chr1  ...  gene_id "DDX11L1"; transcript_id "NR_046018.2"...

>>> start = a.named_gtf.loc["TP63"]["start"].min()
>>> end = a.named_gtf.loc["TP63"]["end"].max()
>>> chrom = a.named_gtf.loc["TP63"]["seqname"][0]

>>> g = genomepy.Genome("hg38", genomes_dir="./data/genomes")
>>> g[chrom][start:end]
>chr3:189596747-189897276
gcaacccgctggggtcaccttccacactgtggaagctttgttcttttgctctttgcagtaaatcttgct...

```

The `genomepy.Genome` class builds on top of the `pyfaidx.Fasta` class, 
see the [pyfaidx documentation](https://github.com/mdshw5/pyfaidx) for more details.
The `genomepy.Annotation` class contains pandas Dataframes with GTF and BED files, as well as additional class methods to utilize these.

## Frequently Asked Questions

Genomepy utilizes external databases to obtain your files.
Unfortunately this sometimes causes issues.
Here are some of the more common issues, with solutions.

Let us know if you encounter issues you cannot solve by creating [a new issue](https://github.com/vanheeringen-lab/genomepy/issues).

### Provider is offline
Occasionally one of the providers experience connection issues, which can last anywhere between minutes to hours.
When this happens genomepy will warn that the provider appears offline, or that the URL seems broken.

If the issue does not pass, you can try to reset genomepy.
Simply run `genomepy clean` on the command line, or run `genomepy.clean()` in Python.

### This genome is missing
Genomepy stores provider data on your computer to rerun it faster later.
If a provider was offline during this time, it may miss (parts of) the data.

To re-download the data, remove the local data with `genomepy clean`, then `search` for your genome again.

### This genome is STILL missing/URL is broken
Sadly, not everything (naming, structure, filenames) is always consistent on the provider end.
Contact the provider to get it fixed!
One notable group are Ensembl fungi, which seems to be mostly mislabelled.

In the meantime, you can still use the power of genomepy by manually retrieving the URLs,
and downloading the files with `genomepy install GENOME_URL -p url --url-to-annotation ANNOTATION_URL`.

### The genomepy config was corrupted
You can create a new one with `genomepy config generate` on command line,
or `genomepy.manage_config("generate")` in Python.

### What's genomepy maximum memory usage?
Genomepy does not read a genome fully into memory. 
Therefore, installing takes less than 1 GB RAM regardless of the genome's size.
Searching NCBI is the most costly operation, using around 3 GB (the first time).

### Which genome/gene annotation to use
Each provider has its pros and cons:
* Ensembl has excellent gene annotations, but their chromosome names can cause issues with some tools.
* UCSC has an excellent genome browser, but their gene annotations vary in format.
* NCBI allows public submissions, and so has the latest versions, although not always complete or error free.

Use `genomepy search` to see your options, and `genomepy annotation` to check the quality of the gene annotation(s).

## Getting help

If you want to report a bug or issue, or have problems with installing or running the software please create [a new issue](https://github.com/vanheeringen-lab/genomepy/issues).
This is the preferred way of getting support.
Alternatively, you can [mail me](mailto:simon.vanheeringen@gmail.com).

## Contributing

Contributions welcome! Send me a pull request or [get in touch](mailto:simon.vanheeringen@gmail.com).

When contributing a PR, please use the [develop](https://github.com/vanheeringen-lab/genomepy/tree/develop) branch.

### Quick development setup: 
1. Fork & download this repo. 
2. `cd` into your local repo. 
3. `git checkout develop`
4. `conda env create -f environment.yaml`
5. `conda activate genomepy`
6. `pip install -e .`
8. `git checkout -b` your_develop_branch

The command line and python imports will now use the code in your local repo. 
To test your changes locally, run the following command: `pytest -vvv`

## Contributors

- Siebren Frölich - [@siebrenf](https://github.com/siebrenf)
- Maarten van der Sande - [@Maarten-vd-Sande](https://github.com/Maarten-vd-Sande)
- Tilman Schäfers [@tilschaef](https://github.com/tilschaef)
- Simon van Heeringen - [@simonvh](https://github.com/simonvh)
- Dohoon Lee - [@dohlee](https://github.com/dohlee)
- Jie Zhu - [@alienzj](https://github.com/alienzj)

## Citation

If you use genomepy in your research, please cite it: [10.21105/joss.00320](http://dx.doi.org/10.21105/joss.00320).

## License

This module is licensed under the terms of the [MIT license](https://opensource.org/licenses/MIT).
