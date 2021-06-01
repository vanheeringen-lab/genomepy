import csv
import os
import re
import subprocess as sp
from tempfile import mkdtemp
from typing import Iterable, Optional, Tuple, Union

from loguru import logger
import mygene
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Memory
from appdirs import user_cache_dir

from genomepy.provider import ProviderBase
from genomepy.utils import (
    get_genomes_dir,
    gzip_and_name,
    gunzip_and_name,
    glob_ext_files,
    read_readme,
    update_readme,
    mkdir_p,
    rm_rf,
    _open,
    best_search_result,
    safe,
)
from genomepy.__about__ import __version__

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
        self.readme_file = os.path.join(self.genome_dir, "README.txt")
        self.annotation_gtf_file = self._get_genome_file("annotation.gtf")
        self.annotation_bed_file = self._get_genome_file("annotation.bed")

        # variables that may be None
        self.genome_file = self._get_genome_file("fa", check_exists=False)
        self.index_file = self.genome_file + ".fai" if self.genome_file else None
        self.sizes_file = self.genome_file + ".sizes" if self.genome_file else None

        # variables loaded on request
        self.genome_contigs = []
        self.annotation_contigs = []
        self.gtf = None
        self.bed = None
        self.genes = []

    def _get_genome_file(self, ext: str, check_exists: Optional[bool] = True):
        """
        Returns the filepath to a single (gzipped) file in the genome_dir with matching ext.
        """
        ext_files = glob_ext_files(self.genome_dir, ext)
        fname = [f for f in ext_files if f"{self.name}.{ext}" in f]
        if len(fname) == 1:
            filepath = os.path.join(self.genome_dir, fname[0])
            return filepath

        # error handling
        if len(fname) > 1:
            raise IndexError(
                f"Too many ({len(fname)}) files match '{self.name}.{ext}(.gz)' "
                f"in directory {self.genome_dir}"
            )
        if check_exists:
            raise FileNotFoundError(
                f"could not find '{self.name}.{ext}(.gz)' in directory {self.genome_dir}"
            )
        return

    @property
    def genome_contigs(self):
        return self.__genome_contigs

    @genome_contigs.setter
    def genome_contigs(self, genome_contig_list):
        self.__genome_contigs = genome_contig_list

    @genome_contigs.getter
    def genome_contigs(self):
        if not self.__genome_contigs:
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
            self.__annotation_contigs = list(set(get_column(self.annotation_bed_file)))
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
            self.__gtf = pd.read_csv(
                self.annotation_gtf_file,
                sep="\t",
                names=GTF_FORMAT,
                dtype={"seqname": str, "start": np.uint32, "end": np.uint32},
                comment="#",
            )
        return self.__gtf

    @property
    def bed(self):
        return self.__bed

    @bed.setter
    def bed(self, bed):
        self.__bed = bed

    @bed.getter
    def bed(self):
        if not isinstance(self.__bed, pd.DataFrame):
            self.__bed = pd.read_csv(
                self.annotation_bed_file,
                sep="\t",
                names=BED12_FORMAT,
                dtype={"chrom": str, "start": np.uint32, "end": np.uint32},
                comment="#",
            )
        return self.__bed

    @property
    def genes(self):
        return self.__genes

    @genes.setter
    def genes(self, genes):
        self.__genes = genes

    @genes.getter
    def genes(self):
        if not self.__genes:
            self.__genes = list(set(self.bed.name))
        return list(set(self.bed.name))

    def gene_coords(self, genes: Iterable[str]) -> pd.DataFrame:
        """
        Retrieve gene locations.

        Parameters
        ----------
        genes : Iterable
            List of gene names as found in the BED file "name" column.

        Returns
        -------
        pandas.DataFrame with gene annotation.
        """
        gene_list = list(genes)
        gene_info = self.bed[["chrom", "start", "end", "name", "strand"]].set_index(
            "name"
        )

        gene_info = gene_info.loc[gene_list]
        if gene_info.shape[0] < 0.9 * len(gene_list):
            logger.warning(
                "Not all genes were found. All gene names can be found in `Annotation.genes`"
            )
        return gene_info.reset_index()[["chrom", "start", "end", "name", "strand"]]

    @staticmethod
    def filter_regex(
        annotation_file: str,
        regex: Optional[str] = ".*",
        invert_match: Optional[bool] = False,
    ):
        """
        Filter annotation file (works with both gtf and bed) by contig name.

        annotation_file: path to annotation file.
        regex: regex string to keep.
        invert_match: keep contigs NOT matching the regex string.

        return a list of discarded contigs.
        """
        old = pd.read_csv(annotation_file, sep="\t", header=None, comment="#")
        pattern = re.compile(regex)
        filter_func = old[0].map(lambda x: bool(pattern.match(x)) is not invert_match)
        new = old[filter_func]
        new.to_csv(
            annotation_file, sep="\t", header=None, index=None, quoting=csv.QUOTE_NONE
        )

        discarded_contigs = list(set(old[~filter_func][0]))
        return discarded_contigs

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

    def gtf_from_bed(self):
        """
        Create an annotation gtf file from an annotation bed file.
        Overwrites the existing gtf file.
        """
        tmp_dir = mkdtemp(dir=self.genome_dir)

        # unzip if needed
        self.annotation_bed_file, is_unzipped = gunzip_and_name(
            self.annotation_bed_file
        )
        new_gtf_file = os.path.join(
            tmp_dir, os.path.basename(re.sub(r"\.gz$", "", self.annotation_gtf_file))
        )

        cmd = "bedToGenePred {0} /dev/stdout | genePredToGtf file /dev/stdin {1}"
        sp.check_call(cmd.format(self.annotation_bed_file, new_gtf_file), shell=True)

        # gzip if needed
        self.annotation_bed_file = gzip_and_name(self.annotation_bed_file, is_unzipped)
        new_gtf_file = gzip_and_name(
            new_gtf_file, self.annotation_gtf_file.endswith(".gz")
        )

        os.replace(new_gtf_file, self.annotation_gtf_file)
        rm_rf(tmp_dir)

    def bed_from_gtf(self):
        """
        Create an annotation bed file from an annotation gtf file.
        Overwrites the existing bed file.
        """
        tmp_dir = mkdtemp(dir=self.genome_dir)

        # unzip if needed
        self.annotation_gtf_file, is_unzipped = gunzip_and_name(
            self.annotation_gtf_file
        )
        new_bed_file = os.path.join(
            tmp_dir, os.path.basename(re.sub(r"\.gz$", "", self.annotation_bed_file))
        )

        cmd = "gtfToGenePred {0} /dev/stdout | genePredToBed /dev/stdin {1}"
        sp.check_call(cmd.format(self.annotation_gtf_file, new_bed_file), shell=True)

        # gzip if needed
        self.annotation_gtf_file = gzip_and_name(self.annotation_gtf_file, is_unzipped)
        new_bed_file = gzip_and_name(
            new_bed_file, self.annotation_bed_file.endswith(".gz")
        )

        os.replace(new_bed_file, self.annotation_bed_file)
        rm_rf(tmp_dir)

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

        match_contigs: attempt to fix contig names?
        filter_contigs: remove contigs from the annotations that are missing from the genome?
        """
        if self.genome_file is None:
            logger.error("A genome is required for sanitizing!")
            return

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
        metadata, _ = read_readme(self.readme_file)
        if metadata["provider"] == "Ensembl":
            return metadata["name"], metadata["assembly_accession"], metadata["tax_id"]

        if metadata["assembly_accession"] == "na":
            logger.warning(
                "Cannot find a matching genome without an assembly accession."
            )
            return

        asm_acc = metadata["assembly_accession"]
        ensembl_search = list(ProviderBase.search_all(asm_acc, provider="Ensembl"))
        search_result = best_search_result(asm_acc, ensembl_search)
        if len(search_result) == 0:
            logger.warning(
                f"No assembly accession found on Ensembl similar to {asm_acc}"
            )
            return

        return safe(search_result[0]), search_result[2], search_result[4]

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

        return query_mygene(query, tax_id, fields, batch_size)

    @staticmethod
    def _filter_query(query: pd.DataFrame) -> pd.DataFrame:
        """per queried gene, keep the best matching, non-NaN, mapping"""
        if "notfound" in query:
            query = query[~query.notfound == np.NaN]  # drop unmatched genes
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
        df["split_id"] = df["name"].str.split(r"\.", expand=True)[0]
        genes = set(df["split_id"])
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

    def _map_locations(self, to: str) -> pd.DataFrame:
        """
        Load chromosome mapping from one version/assembly to another using the
        NCBI assembly reports.

        Parameters
        ----------
        to: str
            target provider (UCSC, Ensembl or NCBI)

        Returns
        -------
        pandas.DataFrame
            Chromosome mapping.
        """
        genomes_dir = os.path.dirname(self.genome_dir)
        mapping = ProviderBase.map_locations(self.name, to, genomes_dir)
        return mapping

    def map_locations(self, annot: Union[str, pd.DataFrame], to: str) -> pd.DataFrame:
        """
        Map chromosome mapping from one version/assembly to another using the
        NCBI assembly reports.

        Drops missing contigs.

        Parameters
        ----------
        annot: pd.Dataframe
            annotation to map (a pandas dataframe, "bed" or "gtf")
        to: str
            target provider (UCSC, Ensembl or NCBI)

        Returns
        -------
        pandas.DataFrame
            Chromosome mapping.
        """
        mapping = self._map_locations(to)
        if isinstance(annot, str):
            df = self.bed if annot == "bed" else self.gtf
        else:
            df = annot
        df = df.set_index(df.columns[0])
        df = df.join(mapping, how="inner")
        df = df.reset_index()
        return df


@memory.cache
def query_mygene(
    query: Iterable[str],
    tax_id: str,
    fields: str = "genomic_pos",
    batch_size: int = 10000,
) -> pd.DataFrame:
    """
    Use mygene.info to map gene identifiers to another type.

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
    batch_size: int, optional
        Controls batch size for REST API.

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
            query[i : i + batch_size],
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
