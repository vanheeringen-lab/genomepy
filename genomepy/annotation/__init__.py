"""Annotation class, modules & related functions"""
import os
import re
from pathlib import Path
from typing import Iterable, Optional, Union

import numpy as np
import pandas as pd
from loguru import logger

from genomepy.annotation.mygene import _map_genes, query_mygene
from genomepy.annotation.sanitize import _sanitize
from genomepy.annotation.utils import _check_property, _parse_annot, read_annot
from genomepy.files import read_readme
from genomepy.providers import map_locations
from genomepy.utils import cleanpath, get_genomes_dir, safe

__all__ = ["Annotation", "query_mygene", "filter_regex"]


class Annotation:
    """
    Manipulate genes and whole gene annotations with pandas dataframes.

    Parameters
    ----------
    name : str
        Genome name/directory/fasta or gene annotation BED/GTF file.
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

    def __init__(self, name: str, genomes_dir: str = None):
        # name and directory
        n, g = _get_name_and_dir(name, genomes_dir)
        self.name = n
        "genome name"
        self.genome_dir = g
        "path to the genome directory"

        # annotation files
        fname = cleanpath(name)
        suffixes = Path(fname).suffixes[-2:]
        b = fname
        if not (".bed" in suffixes or ".BED" in suffixes):
            b = _get_file(self.genome_dir, f"{self.name}.annotation.bed")
        self.annotation_bed_file = b
        "path to the gene annotation BED file"
        g = fname
        if not (".gtf" in suffixes or ".GTF" in suffixes):
            g = _get_file(self.genome_dir, f"{self.name}.annotation.gtf")
        self.annotation_gtf_file = g
        "path to the gene annotation GTF file"

        # genome files
        g = fname
        if ".fa" not in suffixes:
            g = _get_file(self.genome_dir, f"{self.name}.fa", False)
        self.genome_file = g
        "path to the genome fasta"
        self.readme_file = _get_file(self.genome_dir, "README.txt", False)
        "path to the README file"
        self.index_file = _get_file(self.genome_dir, f"{self.name}.fa.fai", False)
        "path to the genome index"
        self.sizes_file = _get_file(self.genome_dir, f"{self.name}.fa.sizes", False)
        "path to the chromosome sizes file"

        # genome attributes
        t = read_readme(str(self.readme_file))[0]["tax_id"]
        self.tax_id = None if t == "na" else int(t)
        "genome taxonomy identifier"

    # lazy attributes
    def __getattribute__(self, name):
        val = super(Annotation, self).__getattribute__(name)
        if val is not None:
            return val

        # if the attribute is None/empty, check if it is a lazy attribute
        if name == "bed":
            _check_property(self.annotation_bed_file, f"{self.name}.annotation.bed")
            val = read_annot(self.annotation_bed_file)
            setattr(self, name, val)

        elif name == "gtf":
            _check_property(self.annotation_gtf_file, f"{self.name}.annotation.gtf")
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
            _check_property(self.sizes_file, f"{self.name}.fa.sizes")
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
        if self.readme_file is None:
            raise AttributeError(
                "Can only map genomepy annotations (a readme file is required)"
            )
        genomes_dir = os.path.dirname(self.genome_dir)
        mapping = map_locations(self.name, to, genomes_dir)
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

    def gtf_dict(
        self, key, value, string_values=True, annot: Union[str, pd.DataFrame] = "gtf"
    ):
        """
        Create a dictionary based on the columns or attribute fields in a GTF.

        Parameters
        ----------
        key : str
            column name or attribute fields (e.g. "seqname", "gene_name")
        value : str
            column name or attribute fields (e.g. "gene_id", "transcript_name")
        string_values : bool, optional
            attempt to format the dict values as strings
            (only happens if all value lists are length 1)
        annot : str or pd.Dataframe, optional
            annotation to filter: "gtf" or a pandas dataframe

        Returns
        -------
        dict
            with values as lists. If string_values is True
            and all lists are length 1, values will be strings.
        """
        df = _parse_annot(self, annot)
        k = key in df.columns
        v = value in df.columns
        # pd.DataFrame.iterrows() is slow. this is not.
        attributes = zip(
            df[key] if k else df.attribute,
            df[value] if v else df.attribute,
        )

        def _get_attr_item(series, item):
            """
            example series.attribute: "...; gene_name: "TP53"; ..."
            item="gene_name" would return "TP53"
            """
            split = series.split(item)
            return split[1].split('"')[1]  # item might not exist

        def _get_col_item(series, _):
            return series

        get_key = _get_col_item if k else _get_attr_item
        get_val = _get_col_item if v else _get_attr_item
        a_dict = dict()
        for row in attributes:
            try:
                k = get_key(row[0], key)
                v = get_val(row[1], value)
            except IndexError:
                continue

            # unique values per key
            if k in a_dict:
                a_dict[k].update({v})
            else:
                a_dict[k] = {v}

        # return values as str if all values are length 1
        # and string_values is True, else return values as list
        all_len_1 = string_values and all(len(v) == 1 for v in a_dict.values())
        for k, v in a_dict.items():
            a_dict[k] = list(v)[0] if all_len_1 else list(v)

        return a_dict


def _get_name_and_dir(name, genomes_dir=None):
    """
    Returns the name and directory of the genome.
    """
    fname = cleanpath(name)
    genomes_dir = get_genomes_dir(genomes_dir, check_exist=False)
    if os.path.isfile(fname):
        exts = ["gtf", "GTF", "bed", "BED", "fa"]
        if not any(ext in fname for ext in exts):
            raise NotImplementedError(
                "Only (gzipped) bed, gtf or fasta files are supported!"
            )
        genome_dir = os.path.dirname(fname)
        name = safe(os.path.basename(fname))
        # remove suffices
        any_ext = "(" + ")|(".join(exts) + ")"
        name = re.sub(fr"(\.annotation)?\.({any_ext})(\.gz)?$", "", name)
    elif os.path.isdir(fname):
        genome_dir = fname
        name = safe(os.path.basename(fname))
    elif name in os.listdir(genomes_dir):
        genome_dir = os.path.join(genomes_dir, name)
    else:
        raise ValueError(f"Could not find {name}")
    return name, genome_dir


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
