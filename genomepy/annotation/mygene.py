import os
from typing import Iterable, Optional, Tuple, Union

import mygene
import pandas as pd
from appdirs import user_cache_dir
from joblib import Memory
from loguru import logger
from tqdm import tqdm

from genomepy.__about__ import __version__
from genomepy.annotation.utils import _check_property
from genomepy.files import read_readme
from genomepy.providers import search_all
from genomepy.utils import best_search_result, mkdir_p, safe

cachedir = os.path.join(user_cache_dir("genomepy"), __version__)
memory = Memory(cachedir, verbose=0)


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


def ensembl_genome_info(self) -> Optional[Tuple[str, str, str]]:
    """
    Return Ensembl genome information for this genome.
    Requires accession numbers to match (excluding patch numbers)

    Returns
    -------
    (str, str, str)
        Ensembl name, accession, taxonomy_id
    """
    _check_property(self.readme_file, "README.txt")

    metadata, _ = read_readme(self.readme_file)
    if metadata.get("provider") == "Ensembl":
        return metadata["name"], metadata["assembly_accession"], metadata["tax_id"]

    if metadata.get("assembly_accession", "na") == "na":
        logger.warning("Cannot find a matching genome without an assembly accession.")
        return

    asm_acc = metadata["assembly_accession"]
    ensembl_search = list(search_all(asm_acc, provider="Ensembl"))
    search_result = best_search_result(asm_acc, ensembl_search)
    if len(search_result) == 0:
        logger.warning(f"No assembly accession found on Ensembl similar to {asm_acc}")
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
    if not str(tax_id).isdigit():
        raise ValueError("No taxomoy ID found")

    return query_mygene(query, tax_id, fields, batch_size)


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
