import csv
import os
import re
from typing import Iterable, Optional, Tuple, Union

import mygene
import numpy as np
import pandas as pd
from appdirs import user_cache_dir
from joblib import Memory
from loguru import logger
from tqdm import tqdm

from genomepy.__about__ import __version__
from genomepy.annotation.utils import _check_property, read_annot
from genomepy.files import _open
from genomepy.providers import map_locations, search_all
from genomepy.utils import best_search_result, get_genomes_dir, mkdir_p, safe


class Annotation:
    """
    Gene annotation related functions

    Parameters
    ----------
    name : str
        Genome name

    genomes_dir : str
        Genomes installation directory

    Returns
    -------
    Annotation object
    """

    from genomepy.annotation.sanitize import sanitize

    def __init__(self, name: str, genomes_dir: str = None):
        self.name = name
        self.genome_dir = os.path.join(get_genomes_dir(genomes_dir), name)

        # annotation files
        self.readme_file = _get_file(self.genome_dir, "README.txt", False)
        self.annotation_gtf_file = _get_file(
            self.genome_dir, f"{self.name}.annotation.gtf"
        )
        self.annotation_bed_file = _get_file(
            self.genome_dir, f"{self.name}.annotation.bed"
        )

        # genome files
        self.genome_file = _get_file(self.genome_dir, f"{self.name}.fa", False)
        self.index_file = _get_file(self.genome_dir, f"{self.name}.fa.fai", False)
        self.sizes_file = _get_file(self.genome_dir, f"{self.name}.fa.sizes", False)

    # lazy attributes
    def __getattr__(self, name):
        if name == "bed":
            _check_property(self.annotation_bed_file, f"{self.name}.annotation.bed")
            val = read_annot(self.annotation_bed_file)
            setattr(self, name, val)

        elif name == "gtf":
            _check_property(self.annotation_gtf_file, f"{self.name}.annotation.gtf")
            val = read_annot(self.annotation_gtf_file)
            setattr(self, name, val)

        elif name == "genome_contigs":
            _check_property(self.sizes_file, f"{self.name}.fa.sizes")
            val = get_column(self.sizes_file)
            setattr(self, name, val)

        elif name == "annotation_contigs":
            val = list(set(self.bed["chrom"]))
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

        else:
            raise AttributeError(
                f"'{self.__class__.__name__}' object has no attribute called '{name}'."
            )

        return getattr(self, name)

    def genes(self, annot: str = "bed") -> list:
        """
        Retrieve gene names from the specified annotation.

        For BED files, names from the 'name' columns

        For GTF files, names from the 'gene_name' field in the attribute column, if available.

        Parameters
        ----------
        annot : str, optional
            Annotation file type: 'bed' or 'gtf'

        Returns
        -------
        list with gene names
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
            List of gene names as found in the given annotation file type.

        annot : str, optional
            Annotation file type: 'bed' or 'gtf'

        Returns
        -------
        pandas.DataFrame with gene annotation.
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
        Map chromosome mapping from one assembly to another using the
        NCBI assembly reports.

        Drops missing contigs.

        Parameters
        ----------
        annot: str or pd.Dataframe
            annotation to map (a pandas dataframe, "bed" or "gtf")
        to: str
            target provider (UCSC, Ensembl or NCBI)
        drop: bool, optional
            if True, replace the chromosome column.
            if False, add a 2nd chromosome column

        Returns
        -------
        pandas.DataFrame
            Chromosome mapping.
        """
        genomes_dir = os.path.dirname(self.genome_dir)
        mapping = map_locations(self.name, to, genomes_dir)
        if mapping is None:
            return

        df = _parse_annot(self, annot)
        if df is None:
            logger.error("Argument 'annot' must be 'gtf', 'bed' or a pandas dataframe.")
            return

        index_name = df.index.name
        if not set([index_name] + df.columns.to_list()) & {"chrom", "seqname"}:
            logger.error(
                "Location mapping requires a column named 'chrom' or 'seqname'."
            )
            return

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
    ) -> Union[pd.DataFrame, None]:
        df = _parse_annot(self, annot)
        if df is None:
            logger.error("Argument 'annot' must be 'gtf', 'bed' or a pandas dataframe.")
            return

        return filter_regex(df, regex, invert_match, column)


def get_column(
    fname: str,
    n: Optional[int] = 0,
    sep: str = "\t",
    comment_char: Union[None, str] = "#",
):
    """
    Return the nth column from a separated file.
    """

    def comment(_: str):
        return False

    if isinstance(comment_char, str):

        def comment(string: str):  # noqa: F811
            return string.startswith(comment_char)

    column = []
    with _open(fname) as f:
        for line in f:
            line = str(line).strip()
            if not comment(line):
                column.append(line.split(sep)[n])
    return column


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


def _parse_annot(self, annot):
    if isinstance(annot, pd.DataFrame):
        df = annot
    elif isinstance(annot, str) and annot == "bed":
        df = self.bed
    elif isinstance(annot, str) and annot == "gtf":
        df = self.gtf
    else:
        df = None
    return df


def filter_regex(
    df: pd.DataFrame,
    regex: str,
    invert_match: Optional[bool] = False,
    column: Union[str, int] = 0,
) -> Union[pd.DataFrame, None]:
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
    filtered pd.DataFrame
    """
    if column not in df.columns:
        if isinstance(column, int):
            column = df.columns[column]
        else:
            logger.error(
                f"Column '{column}' not found in annotation columns {list(df.columns)}"
            )
            return

    pattern = re.compile(regex)
    filter_func = df[column].map(lambda x: bool(pattern.match(x)) is not invert_match)
    return df[filter_func]
