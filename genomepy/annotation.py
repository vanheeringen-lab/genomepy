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
from genomepy.files import _open, generate_annot, read_readme, update_readme
from genomepy.provider import Provider
from genomepy.utils import best_search_result, get_genomes_dir, mkdir_p, safe

GTF_FORMAT = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]

BED12_FORMAT = [
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRrgb",
    "blockCount",
    "blockSizes",
    "blockStarts",
]

cachedir = os.path.join(user_cache_dir("genomepy"), __version__)
memory = Memory(cachedir, verbose=0)


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

    def __init__(self, name: str, genomes_dir: str = None):
        self.name = name
        self.genome_dir = os.path.join(get_genomes_dir(genomes_dir), name)
        self.readme_file = self._get_file("README.txt", warn_missing=False)
        self.annotation_gtf_file = self._get_file(f"{self.name}.annotation.gtf")
        self.annotation_bed_file = self._get_file(f"{self.name}.annotation.bed")

        # genome files
        self.genome_file = self._get_file(f"{self.name}.fa", warn_missing=False)
        self.index_file = self._get_file(f"{self.name}.fa.fai", warn_missing=False)
        self.sizes_file = self._get_file(f"{self.name}.fa.sizes", warn_missing=False)

        # variables loaded on request
        self.genome_contigs = []
        self.annotation_contigs = []
        self.gtf = None
        self.named_gtf = None  # HGNC names on index
        self.bed = None

    def _get_file(self, fname: str, warn_missing: Optional[bool] = True):
        """
        Returns the filepath to a single (gzipped) file in the genome_dir with matching ext.
        """
        fpath = os.path.join(self.genome_dir, fname)
        if os.path.exists(fpath):
            return fpath
        if os.path.exists(f"{fpath}.gz"):
            return f"{fpath}.gz"
        if warn_missing:
            logger.warning(
                f"Could not find '{fname}(.gz)' in directory {self.genome_dir}. "
                "Methods using this file won't work!"
            )
        return

    @staticmethod
    def _check_required(prop, fname):
        if prop is None:
            raise FileNotFoundError(f"'{fname}' required.")

    @property
    def genome_contigs(self):
        return self.__genome_contigs

    @genome_contigs.setter
    def genome_contigs(self, genome_contig_list):
        self.__genome_contigs = genome_contig_list

    @genome_contigs.getter
    def genome_contigs(self):
        if not self.__genome_contigs:
            self._check_required(self.sizes_file, f"{self.name}.fa.sizes")
            self.__genome_contigs = get_column(self.sizes_file)
        return self.__genome_contigs

    @property
    def annotation_contigs(self):
        return self.__annotation_contigs

    @annotation_contigs.setter
    def annotation_contigs(self, annotation_contig_list):
        self.__annotation_contigs = annotation_contig_list

    @annotation_contigs.getter
    def annotation_contigs(self):
        if not self.__annotation_contigs:
            self.__annotation_contigs = list(set(self.bed["chrom"]))
        return self.__annotation_contigs

    @property
    def gtf(self):
        return self.__gtf

    @gtf.setter
    def gtf(self, gtf):
        self.__gtf = gtf

    @gtf.getter
    def gtf(self):
        if not isinstance(self.__gtf, pd.DataFrame):
            self._check_required(
                self.annotation_gtf_file, f"{self.name}.annotation.gtf"
            )
            self.__gtf = pd.read_csv(
                self.annotation_gtf_file,
                sep="\t",
                names=GTF_FORMAT,
                dtype={"seqname": str, "start": np.uint32, "end": np.uint32},
                comment="#",
            )
        return self.__gtf

    @property
    def named_gtf(self):
        return self.__named_gtf

    @named_gtf.setter
    def named_gtf(self, gtf):
        self.__named_gtf = gtf

    @named_gtf.getter
    def named_gtf(self):
        """
        GTF dataframe with attribute "gene_name" as index.

        Drops rows without gene_name in the attribute field.
        """
        if not isinstance(self.__named_gtf, pd.DataFrame):
            df = self.gtf[self.gtf.attribute.str.contains("gene_name")]
            names = []
            for row in df.attribute:
                name = str(row).split("gene_name")[1].split(";")[0]
                names.append(name.replace('"', "").replace(" ", ""))
            df["gene_name"] = names
            self.__named_gtf = df.set_index("gene_name")
        return self.__named_gtf

    @property
    def bed(self):
        return self.__bed

    @bed.setter
    def bed(self, bed):
        self.__bed = bed

    @bed.getter
    def bed(self):
        if not isinstance(self.__bed, pd.DataFrame):
            self._check_required(
                self.annotation_bed_file, f"{self.name}.annotation.bed"
            )
            self.__bed = pd.read_csv(
                self.annotation_bed_file,
                sep="\t",
                names=BED12_FORMAT,
                dtype={"chrom": str, "start": np.uint32, "end": np.uint32},
                comment="#",
            )
        return self.__bed

    def genes(self, annot: str = "bed") -> list:
        """
        Retrieve gene names from the specified annotation.

        For BED files, the output names vary, but is always available.

        For GTF files, the output is always HGNC names, if available.

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
            List of gene names as found in the BED file "name" column.

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
        self,
        annot: Union[str, pd.DataFrame],
        regex: Optional[str] = ".*",
        invert_match: Optional[bool] = False,
        column: Union[str, int] = 0,
    ) -> Union[pd.DataFrame, None]:
        """
        Filter a pandas dataframe by a column (default: 1st, contig name).

        Parameters
        ----------
        annot: str or pd.Dataframe
            annotation to map (a pandas dataframe, "bed" or "gtf")

        regex : str
            regex string to keep

        invert_match : bool, optional
            keep contigs NOT matching the regex string

        column: str or int, optional
            column name or number to filter (default: 1st, contig name)

        Returns
        -------
        filtered pandas.DataFrame
        """
        df = self._parse_annot(annot)
        if df is None:
            logger.error("Argument 'annot' must be 'gtf', 'bed' or a pandas dataframe.")
            return

        if column not in df.columns:
            if isinstance(column, int):
                column = df.columns[column]
            else:
                logger.error(
                    f"Column '{column}' not found in annotation columns {list(df.columns)}"
                )
                return

        pattern = re.compile(regex)
        filter_func = df[column].map(
            lambda x: bool(pattern.match(x)) is not invert_match
        )
        return df[filter_func]

    def _filter_genome_contigs(self):
        """
        filter the gtf for contigs in the genome.
        return a list of discarded contigs.
        """
        old = pd.read_csv(self.annotation_gtf_file, sep="\t", header=None, comment="#")
        filter_func = old[0].isin(self.genome_contigs)
        new = old[filter_func]
        new.to_csv(
            self.annotation_gtf_file,
            sep="\t",
            header=None,
            index=None,
            quoting=csv.QUOTE_NONE,
        )

        discarded_contigs = list(set(old[~filter_func][0]))
        return discarded_contigs

    def _is_conforming(self):
        """
        Check if genome and annotation contigs conform.
        Returns bool.
        """
        overlapping_contigs = set(self.genome_contigs) & set(self.annotation_contigs)
        return bool(overlapping_contigs)

    def _conforming_index(self) -> int:
        """
        Check if the annotation contigs can be made to conform the genome contigs.
        Returns int: -1 = not conformable, 0 = conforming, >0 = index of conforming element.
        """
        # example conformable files
        # genome.fa:
        #
        # >chr1 alternate_name_chr1 more_info                    <-- header
        # ATCGATCGATCGATCG                                       <-- sequence
        #
        # genome.annotation.gtf:
        #
        # seqname              source    feature start  end ...  <-- header
        # alternate_name_chr1  provider  exon    0      100 ...
        #
        # returns 1 because genome header[1] matches the GTF contig names

        # all headers in the file should be formatted the same. so take the first.
        header = []
        with _open(self.genome_file, "r") as fa:
            for line in fa:
                if line.startswith(">"):
                    header = line.strip(">\n").split(" ")
                    break

        for contig in self.annotation_contigs:
            if contig in header:
                matching_contig_index = header.index(contig)
                return matching_contig_index
        return -1

    @staticmethod
    def _is_conformable(matching_contig_index: int):
        if matching_contig_index == 0:
            raise ValueError("genome and annotation are conforming.")
        return matching_contig_index > 0

    def _contig_conversion_dict(self, matching_element: int):
        """
        build a conversion dict
        returns the dict and a list of duplicate contig names (should be empty!)
        """
        conversion_table = {}
        duplicate_contigs = []
        with _open(self.genome_file, "r") as fa:
            for line in fa:
                if line.startswith(">"):
                    line = line.strip(">\n").split(" ")
                    if line[matching_element] in conversion_table:
                        duplicate_contigs.append(line[matching_element])
                    conversion_table[line[matching_element]] = line[0]
        return conversion_table, list(set(duplicate_contigs))

    def _conform_gtf(
        self, conversion_dict: dict, filter_contigs: Optional[bool] = True
    ):
        """
        Create an annotation gtf file with contigs conforming to the genome
        Overwrites the existing gtf file.
        Returns a list of contigs not present in the genome
        """
        old = pd.read_csv(self.annotation_gtf_file, sep="\t", header=None, comment="#")
        converted_contigs = old[0].map(conversion_dict)  # contigs missing become NaNs
        nans = converted_contigs.isna()

        new = old.copy()
        new[0] = converted_contigs.fillna(old[0])  # convert contigs (keeps missing)
        if filter_contigs:
            new = new[~nans]  # remove missing contigs
        new.to_csv(
            self.annotation_gtf_file,
            sep="\t",
            header=None,
            index=None,
            quoting=csv.QUOTE_NONE,
        )

        missing_contigs = list(set(old[nans][0]))
        return missing_contigs

    def bed_from_gtf(self):
        """
        Create an annotation bed file from an annotation gtf file.
        Overwrites the existing bed file.
        """
        generate_annot(
            self.annotation_gtf_file, self.annotation_bed_file, overwrite=True
        )

    def gtf_from_bed(self):
        """
        Create an annotation gtf file from an annotation bed file.
        Overwrites the existing gtf file.
        """
        generate_annot(
            self.annotation_bed_file, self.annotation_gtf_file, overwrite=True
        )

    def sanitize(
        self,
        match_contigs: Optional[bool] = True,
        filter_contigs: Optional[bool] = True,
    ):
        """
        Compares the genome and gene annotations.
        If the contig names do not conform, attempt to fix this.
        If the contig names conform (after fixing), filter the annotation contigs for genome contigs.
        Log the results in the README.

        Parameters
        ----------
        match_contigs: bool, optional
            attempt to fix contig names?

        filter_contigs: bool, optional
            remove contigs from the annotations that are missing from the genome?
        """
        if self.genome_file is None:
            logger.error("A genome is required for sanitizing!")
            return
        self._check_required(self.readme_file, "README.txt")

        extra_lines = []
        if self._is_conforming():
            status = "contigs match but not filtered"
            if filter_contigs:
                status = "contigs match"
                contigs_filtered_out = self._filter_genome_contigs()
                if contigs_filtered_out:
                    self.bed_from_gtf()
                    status = "contigs match and filtered"
                    extra_lines = [
                        "",
                        "The following contigs were filtered out of the gene annotation:",
                        f"{', '.join(contigs_filtered_out)}",
                    ]
            update_readme(
                self.readme_file, {"sanitized annotation": status}, extra_lines
            )
            return

        if match_contigs is False:
            return

        matching_contig_index = self._conforming_index()
        if self._is_conformable(matching_contig_index) is False:
            # example not possible: ASM2732v1
            update_readme(self.readme_file, {"sanitized annotation": "not possible"})
            return

        conversion_dict, duplicate_contigs = self._contig_conversion_dict(
            matching_contig_index
        )
        if duplicate_contigs:
            logger.warning(
                "\nThe genome contains duplicate contig names, this should really not happen!\n"
                "You may wish to consider filtering these out, or look for another genome assembly..."
                "The following duplicate contigs were found: "
                f"{', '.join(duplicate_contigs)}.\n"
            )
        missing_contigs = self._conform_gtf(conversion_dict, filter_contigs)
        self.bed_from_gtf()
        status = "contigs fixed"
        if missing_contigs:
            status = "contigs fixed and filtered"
            extra_lines = [
                "",
                "The following contigs were filtered out of the gene annotation:",
                f"{', '.join(missing_contigs)}",
            ]
            if filter_contigs is False:
                status = "contigs fixed but not filtered"
                extra_lines[1] = (
                    "The following contigs could not be sanitized"
                    " in the gene annotation, and were kept as-is:"
                )
            logger.info(" ".join(extra_lines))
        update_readme(self.readme_file, {"sanitized annotation": status}, extra_lines)

    def ensembl_genome_info(self) -> Optional[Tuple[str, str, str]]:
        """
        Return Ensembl genome information for this genome.
        Requires accession numbers to match (excluding patch numbers)

        Returns
        -------
        (str, str, str)
            Ensembl name, accession, taxonomy_id
        """
        self._check_required(self.readme_file, "README.txt")

        metadata, _ = read_readme(self.readme_file)
        if metadata.get("provider") == "Ensembl":
            return metadata["name"], metadata["assembly_accession"], metadata["tax_id"]

        if metadata.get("assembly_accession", "na") == "na":
            logger.warning(
                "Cannot find a matching genome without an assembly accession."
            )
            return

        asm_acc = metadata["assembly_accession"]
        ensembl_search = list(Provider.search_all(asm_acc, provider="Ensembl"))
        search_result = best_search_result(asm_acc, ensembl_search)
        if len(search_result) == 0:
            logger.warning(
                f"No assembly accession found on Ensembl similar to {asm_acc}"
            )
            return

        return safe(search_result[0]), search_result[2], search_result[3]

    def _query_mygene(
        self,
        query: Iterable[str],
        fields: str = "genomic_pos",
        batch_size: int = 5000,
    ) -> pd.DataFrame:
        # mygene.info only queries the most recent version of the Ensembl database
        # We can only safely continue if the local genome matched the Ensembl genome.
        # Even if the local genome was installed via Ensembl, we still need to check
        # if it is the same version
        ensembl_info = self.ensembl_genome_info()
        if ensembl_info is None:
            return pd.DataFrame()
        _, _, tax_id = ensembl_info
        if not isinstance(tax_id, int):
            raise ValueError("No taxomoy ID found")

        return query_mygene(query, tax_id, fields, batch_size)

    @staticmethod
    def _filter_query(query: pd.DataFrame) -> pd.DataFrame:
        """per queried gene, keep the best matching, non-NaN, mapping"""
        if "notfound" in query:
            query = query[query.notfound.astype(str) != "nan"]  # drop unmatched genes
            query = query.drop(columns="notfound")
        query = query.dropna()
        if "query" in query:
            query = query.groupby("query").first()  # already sorted by mapping score
            query = query.drop(columns=["_id", "_score"])
        return query

    def map_with_mygene(self, query: Iterable[str], fields: str = "genomic_pos"):
        """
        Use mygene.info to map gene identifiers to another type.

        If the identifier can't be mapped, it will be dropped from the resulting
        annotation. If multiple identifiers match, the best match is used.

        Parameters
        ----------
        query: iterable
            a list or list-like of gene identifiers
        fields : str, optional
            Target identifier to map the query genes to. Valid fields
            are: ensembl.gene, entrezgene, symbol, name, refseq, entrezgene. Note that
            refseq will return the protein refseq_id by default, use `product="rna"` to
            return the RNA refseq_id. Currently, mapping to Ensembl transcript ids is
            not supported.

        Returns
        -------
        pandas.DataFrame with mapped gene annotation.
        """
        # load mapping
        gene_hash = hash(tuple(set(query)))  # hash unique names
        mygene_mapping = os.path.join(
            self.genome_dir, "mygene", f"{fields}_genehash_{gene_hash}.tsv"
        )
        if os.path.exists(mygene_mapping):
            return pd.read_csv(mygene_mapping, sep="\t", index_col=0)

        # request mapping
        ret = self._query_mygene(query, fields)
        ret = self._filter_query(ret)

        # save mapping
        mkdir_p(os.path.join(self.genome_dir, "mygene"))
        ret.to_csv(mygene_mapping, sep="\t", index=True, header=True)
        return ret

    def _map_genes(
        self, df: pd.DataFrame, to: str, product: str = "protein"
    ) -> pd.DataFrame:
        cols = df.columns  # starting columns
        # remove version numbers from gene IDs
        split_id = df["name"].str.split(r"\.", expand=True)[0]
        df = df.assign(split_id=split_id.values)
        genes = set(split_id)
        result = self.map_with_mygene(genes, fields=to)
        if len(result) == 0:
            logger.error("Could not map using mygene.info")
            return pd.DataFrame()
        df = df.join(result, on="split_id")

        # Only in case of RefSeq annotation the product needs to be specified.
        if to == "refseq":
            to = f"{to}.translation.{product}"

        # Get rid of extra columns from query
        df["name"] = df[to]
        df = df[cols].dropna()
        return df

    def map_genes(
        self, gene_field: str, product: str = "protein", df: pd.DataFrame = None
    ) -> Optional[pd.DataFrame]:
        """
        Uses mygene.info to map gene identifiers on the fly by specifying
        `gene_field`. If the identifier can't be mapped, it will be dropped
        from the resulting annotation.

        Returns the dataframe (BED12 by default) with remapped "name" column.

        Parameters
        ----------
        gene_field : str, optional
            Identifier for gene annotation. Uses mygene.info to map ids. Valid fields
            are: ensembl.gene, entrezgene, symbol, name, refseq, entrezgene. Note that
            refseq will return the protein refseq_id by default, use `product="rna"` to
            return the RNA refseq_id. Currently, mapping to Ensembl transcript ids is
            not supported.
        product : str, optional
            Either "protein" or "rna". Only used when `gene_field="refseq"`
        df : pandas.Dataframe, optional
            dataframe with "name" column to remap, will use the annotation's BED file if unspecified.

        Returns
        -------
        pandas.DataFrame with gene annotation.
        """
        if df is None:
            df = self.bed.copy()

        product = product.lower()
        if product not in ["rna", "protein"]:
            raise ValueError("Argument product should be either 'rna' or 'protein'")

        gene_field = gene_field.lower()
        if gene_field not in [
            "ensembl.gene",
            "entrezgene",
            "symbol",
            "name",
            "refseq",
            "entrezgene",
        ]:
            raise ValueError(
                "Argument product should be either 'ensembl.gene', "
                "'entrezgene,' 'symbol', 'name', 'refseq' or 'entrezgene'"
            )

        df = self._map_genes(df, gene_field, product)
        return df

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
        mapping = Provider.map_locations(self.name, to, genomes_dir)
        if mapping is None:
            return

        df = self._parse_annot(annot)
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


@memory.cache
def query_mygene(
    query: Iterable[str],
    tax_id: Union[str, int],
    fields: str = "genomic_pos",
    batch_size: int = 10000,
) -> pd.DataFrame:
    """
    Use mygene.info to map gene identifiers to another type.

    Parameters
    ----------
    query: iterable
        a list or list-like of gene identifiers
    tax_id: str or int
        Target genome taxonomy id
    fields : str, optional
        Target identifier to map the query genes to. Valid fields
        are: ensembl.gene, entrezgene, symbol, name, refseq, entrezgene. Note that
        refseq will return the protein refseq_id by default, use `product="rna"` to
        return the RNA refseq_id. Currently, mapping to Ensembl transcript ids is
        not supported.
    batch_size: int, optional
        Controls batch size for the REST API.

    Returns
    -------
    pandas.DataFrame with mapped gene annotation.
    """
    logger.info("Querying mygene.info...")
    query = [_ for _ in query]
    query_len = len(query)
    it = range(0, query_len, batch_size)
    if query_len > batch_size:
        logger.info("Large query, running in batches...")
        it = tqdm(it)

    result = pd.DataFrame()
    for i in it:
        mg = mygene.MyGeneInfo()
        _result = mg.querymany(
            query[i : i + batch_size],  # noqa
            scopes="symbol,name,ensembl.gene,entrezgene,ensembl.transcript,ensembl,accession.protein,accession.rna",
            fields=fields,
            species=tax_id,
            as_dataframe=True,
            verbose=False,
        )
        result = pd.concat((result, _result))

    if "notfound" in result and result.shape[1] == 1:
        logger.warning("No matching genes found")

    return result


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
