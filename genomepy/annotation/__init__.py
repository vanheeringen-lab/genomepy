"""Annotation class, modules & related functions"""
import os
import re
from pathlib import Path
from typing import Iterable, Optional, Union

import numpy as np
import pandas as pd
from loguru import logger

from genomepy.annotation.mygene import map_genes as _map_genes
from genomepy.annotation.mygene import query_mygene
from genomepy.annotation.sanitize import sanitize as _sanitize
from genomepy.annotation.utils import _check_property, _parse_annot, read_annot
from genomepy.providers import map_locations
from genomepy.utils import get_genomes_dir

__all__ = ["Annotation", "query_mygene", "filter_regex"]


class Annotation:
    """
    Manipulate genes and whole gene annotations with pandas dataframes.

    Parameters
    ----------
    genome : str
        Genome name.
    name : str, optional
        Name of annotation file.
        If name is not specified the default annotation for the genome is used.
    genomes_dir : str, optional
        Genomes installation directory.

    Returns
    -------
    object
        attributes & methods to manipulate gene annotations
    """

    # import methods
    map_genes = _map_genes
    sanitize = _sanitize

    # lazy attributes (loaded when called)
    # listed here for code autocompletion
    bed: pd.DataFrame = None
    "Dataframe with BED format annotation"
    gtf: pd.DataFrame = None
    "Dataframe with GTF format annotation"
    named_gtf: pd.DataFrame = None
    "Dataframe with GTF format annotation, with gene_name as index"
    genome_contigs: list = None
    "Contigs found in the genome fasta"
    annotation_contigs: list = None
    "Contigs found in the gene annotation BED"

    def __init__(self, genome: str, name: str = None, genomes_dir: str = None):
        self.genome = genome
        self.genome_dir = os.path.join(get_genomes_dir(genomes_dir), genome)
        if not os.path.exists(self.genome_dir):
            raise ValueError(f"Genome {self.genome} not found!")

        # annotation file provided
        if name:
            suffixes = Path(name).suffixes[-2:]
            if ".bed" in suffixes or ".BED" in suffixes:
                self.annotation_bed_file = name
            elif ".gtf" in suffixes or ".GTF" in suffixes:
                self.annotation_gtf_file = name
            else:
                raise NotImplementedError(
                    "Only (gzipped) bed and gtf files are supported at the moment!"
                )
        else:
            # annotation files
            self.annotation_gtf_file = _get_file(
                self.genome_dir, f"{self.genome}.annotation.gtf"
            )
            self.annotation_bed_file = _get_file(
                self.genome_dir, f"{self.genome}.annotation.bed"
            )

        # genome files
        self.readme_file = _get_file(self.genome_dir, "README.txt", False)
        self.genome_file = _get_file(self.genome_dir, f"{self.genome}.fa", False)
        self.index_file = _get_file(self.genome_dir, f"{self.genome}.fa.fai", False)
        self.sizes_file = _get_file(self.genome_dir, f"{self.genome}.fa.sizes", False)

    # lazy attributes
    def __getattribute__(self, name):
        val = super(Annotation, self).__getattribute__(name)
        if val is not None:
            return val

        # if the attribute is None/empty, check if it is a lazy attribute
        if name == "bed":
            _check_property(self.annotation_bed_file, f"{self.genome}.annotation.bed")
            val = read_annot(self.annotation_bed_file)
            setattr(self, name, val)

        elif name == "gtf":
            _check_property(self.annotation_gtf_file, f"{self.genome}.annotation.gtf")
            val = read_annot(self.annotation_gtf_file)
            setattr(self, name, val)

        elif name == "named_gtf":
            df = self.gtf[self.gtf.attribute.str.contains("gene_name")]
            names = []
            for row in df.attribute:
                name = str(row).split("gene_name")[1].split(";")[0]
                names.append(name.replace('"', "").replace(" ", ""))
            df = df.assign(gene_name=names)
            val = df.set_index("gene_name")
            setattr(self, name, val)

        elif name == "genome_contigs":
            _check_property(self.sizes_file, f"{self.genome}.fa.sizes")
            val = list(
                set(pd.read_csv(self.sizes_file, sep="\t", header=None, dtype=str)[0])
            )
            setattr(self, name, val)

        elif name == "annotation_contigs":
            val = list(set(self.bed.chrom))
            setattr(self, name, val)

        return val

    # lazily update attributes if upstream attribute is updated
    def __setattr__(self, name, value):
        if name == "bed":
            self.annotation_contigs = None  # noqa
        elif name == "gtf":
            self.named_gtf = None  # noqa
        elif name == "sizes_file":
            self.genome_contigs = None  # noqa
        super(Annotation, self).__setattr__(name, value)

    def genes(self, annot: str = "bed") -> list:
        """
        Retrieve gene names from an annotation.

        For BED files, names are taken from the 'name' columns.

        For GTF files, names are taken from the 'gene_name' field
        in the attribute column, if available.

        Parameters
        ----------
        annot : str, optional
            Annotation file type: 'bed' or 'gtf' (default: "bed")

        Returns
        -------
        list
            gene names
        """
        if annot.lower() == "bed":
            return list(set(self.bed.name))
        return list(set(self.named_gtf.index))

    def gene_coords(self, genes: Iterable[str], annot: str = "bed") -> pd.DataFrame:
        """
        Retrieve gene locations.

        Parameters
        ----------
        genes : Iterable
            List of gene names as found in the given annotation file type
        annot : str, optional
            Annotation file type: 'bed' or 'gtf' (default: "bed")

        Returns
        -------
        pandas.DataFrame
            gene annotation
        """
        gene_list = list(genes)
        if annot.lower() == "bed":
            df = self.bed.set_index("name")
            gene_info = df[["chrom", "start", "end", "strand"]]
        else:
            df = self.named_gtf
            # 1 row per gene
            df = (
                df.groupby(["gene_name", "seqname", "strand"])
                .agg({"start": np.min, "end": np.max})
                .reset_index(level=["seqname", "strand"])
            )
            gene_info = df[["seqname", "start", "end", "strand"]]

        gene_info = gene_info.reindex(gene_list).dropna()
        pct = int(100 * len(set(gene_info.index)) / len(gene_list))
        if pct < 90:
            logger.warning(
                (f"Only {pct}% of genes was found. " if pct else "No genes found. ")
                + "A list of all gene names can be found with `Annotation.genes()`"
            )

        if annot.lower() == "bed":
            return gene_info.reset_index()[["chrom", "start", "end", "name", "strand"]]
        else:
            return gene_info.reset_index()[
                ["seqname", "start", "end", "gene_name", "strand"]
            ]

    def map_locations(
        self, annot: Union[str, pd.DataFrame], to: str, drop=True
    ) -> Union[None, pd.DataFrame]:
        """
        Map chromosome mapping from one assembly to another.

        Uses the NCBI assembly reports to find contigs.
        Drops missing contigs.

        Parameters
        ----------
        annot : str or pd.Dataframe
            annotation to map: "bed", "gtf" or a pandas dataframe.
        to: str
            target provider (UCSC, Ensembl or NCBI)
        drop: bool, optional
            if True, replace the chromosome column.
            If False, add a 2nd chromosome column.

        Returns
        -------
        pandas.DataFrame
            chromosome mapping.
        """
        genomes_dir = os.path.dirname(self.genome_dir)
        mapping = map_locations(self.genome, to, genomes_dir)
        if mapping is None:
            return

        df = _parse_annot(self, annot)
        index_name = df.index.name
        if not set([index_name] + df.columns.to_list()) & {"chrom", "seqname"}:
            raise ValueError(
                "Location mapping requires a column named 'chrom' or 'seqname'."
            )

        # join mapping on chromosome column and return with original index
        is_indexed = df.index.to_list() != list(range(df.shape[0]))
        if is_indexed:
            df = df.reset_index(level=index_name)
        index_col = "chrom" if "chrom" in df.columns else "seqname"
        df = df.set_index(index_col)
        df = mapping.join(df, how="inner")
        df = df.reset_index(drop=drop)
        df.columns = [index_col] + df.columns.to_list()[1:]
        if is_indexed:
            df = df.set_index(index_name if index_name else "index")

        return df

    def filter_regex(
        self,
        annot: Union[str, pd.DataFrame],
        regex: Optional[str] = ".*",
        invert_match: Optional[bool] = False,
        column: Union[str, int] = 0,
    ) -> pd.DataFrame:
        """
        Filter a dataframe by any column using regex.

        Parameters
        ----------
        annot : str or pd.Dataframe
            annotation to filter: "bed", "gtf" or a pandas dataframe
        regex : str
            regex string to match
        invert_match : bool, optional
            keep contigs NOT matching the regex string
        column: str or int, optional
            column name or number to filter (default: 1st, contig name)

        Returns
        -------
        pd.DataFrame
            filtered dataframe
        """
        df = _parse_annot(self, annot)
        return filter_regex(df, regex, invert_match, column)


def _get_file(genome_dir: str, fname: str, warn_missing: Optional[bool] = True):
    """
    Returns the filepath to a single (gzipped) file in the genome_dir with matching ext.
    """
    fpath = os.path.join(genome_dir, fname)
    if os.path.exists(fpath):
        return fpath
    if os.path.exists(f"{fpath}.gz"):
        return f"{fpath}.gz"
    if warn_missing:
        logger.warning(
            f"Could not find '{fname}(.gz)' in directory {genome_dir}. "
            "Methods using this file won't work!"
        )
    return


def filter_regex(
    df: pd.DataFrame,
    regex: str,
    invert_match: Optional[bool] = False,
    column: Union[str, int] = 0,
) -> pd.DataFrame:
    """
    Filter a pandas dataframe by a column (default: 1st, contig name).

    Parameters
    ----------
    df: pd.Dataframe
        annotation to filter (a pandas dataframe)
    regex : str
        regex string to match
    invert_match : bool, optional
        keep contigs NOT matching the regex string
    column: str or int, optional
        column name or number to filter (default: 1st, contig name)

    Returns
    -------
    pd.DataFrame
        filtered dataframe
    """
    if column not in df.columns:
        if isinstance(column, int):
            column = df.columns[column]
        else:
            raise ValueError(
                f"Column '{column}' not found in annotation columns {list(df.columns)}"
            )

    pattern = re.compile(regex)
    filter_func = df[column].map(lambda x: bool(pattern.match(x)) is not invert_match)
    return df[filter_func]
