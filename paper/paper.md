---
title: 'genomepy: download genomes the easy way'
tags:
  - genomics
  - bioinformatics
  - computational biology
authors:
 - name: Simon J. van Heeringen
   orcid: 0000-0002-0411-3219
   affiliation: 1
affiliations:
 - name: Radboud University, Nijmegen, the Netherlands
   index: 1
date: 10 July 2017
bibliography: paper.bib
---

# Summary

The technical advances in high-throughput DNA sequencing have resulted in data
sets of increasing size and complexity. It is essential that experiments are
analyzed in a reproducible manner. Computational workflow systems allow for
automation of analysis pipelines, that scale from personal computers to HPC and
cloud environments [@Koster2012a;@DiTommaso2017;@Vivian2017]. While
standardized solutions exist for installation of software and analysis tools
[@pip;@bioconda], data download often has to either be performed manually or
be scripted on a case-by-case basis.

In many analysis workflows, sequencing reads will need to be mapped to a
reference genome. Here we present genomepy, a simple software package to
automate the download of genomic sequences. It contains both command-line
tools as well as a Python API. Supported providers for genomes include UCSC,
NCBI and Ensembl. Downloaded genome sequences can be soft- or hard-masked and
specific chromosomes or scaffolds can be either included or excluded based on
regular expressions. Genomepy is free and open source software and can be
installed through standard package managers [@pip;@bioconda].

In short, genomepy enables simple, straightforward and automatic downloads of
genomic sequences.

# References
