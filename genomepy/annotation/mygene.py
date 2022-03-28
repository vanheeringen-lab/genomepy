"""Mygene.info functions and methods"""
from typing import Iterable, Union

import mygene
import pandas as pd
from loguru import logger
from tqdm import tqdm

from genomepy.annotation.utils import _parse_annot
from genomepy.caching import memory


def _map_genes(
    self,  # noqa
    field: str,
    product: str = "protein",
    annot: Union[str, pd.DataFrame] = "bed",
) -> pd.DataFrame:
    """
    Use mygene.info to map gene identifiers to any specified `field`.

    Returns the dataframe with remapped "name" column.
    Drops missing identifiers.

    Parameters
    ----------
    annot: str or pd.Dataframe
        Annotation dataframe to map (a pandas dataframe or "bed").
        Is mapped to a column named "name" (required).
    field : str, optional
        Identifier for gene annotation. Uses mygene.info to map ids. Valid fields
        are: ensembl.gene, entrezgene, symbol, name, refseq, entrezgene. Note that
        refseq will return the protein refseq_id by default, use `product="rna"` to
        return the RNA refseq_id. Currently, mapping to Ensembl transcript ids is
        not supported.
    product : str, optional
        Either "protein" or "rna". Only used when `field="refseq"`

    Returns
    -------
    pandas.DataFrame
        remapped gene annotation
    """
    if not isinstance(self.tax_id, int):
        raise AttributeError(
            "A taxonomy identifier is required. "
            "You can set 'Annotation.tax_id' manually"
        )
    to, product = _parse_mygene_input(field, product)
    df = _parse_annot(self, annot)
    if df is None:
        raise ValueError("Argument 'annot' must be 'bed' or a pandas dataframe.")

    cols = df.columns  # starting columns
    if "name" not in cols:
        raise ValueError("Column 'name' is required to map to.")

    # remove version numbers from gene IDs
    split_id = df["name"].str.split(r"\.", expand=True)[0]
    df = df.assign(split_id=split_id.values)
    genes = set(split_id)

    result = query_mygene(genes, self.tax_id, field)
    # result = _query_mygene(self, genes, field=to)
    result = _filter_query(result)
    if len(result) == 0:
        logger.warning("Could not map using mygene.info")
        return pd.DataFrame()
    df = df.join(result, on="split_id")

    # Only in case of RefSeq annotation the product needs to be specified.
    if to == "refseq":
        to = f"{to}.translation.{product}"

    # Get rid of extra columns from query
    df = df.assign(name=df[to].values)  # df["name"] = df[to]
    df = df[cols].dropna()
    return df


@memory.cache
def query_mygene(
    query: Iterable[str],
    tax_id: Union[str, int],
    field: str = "genomic_pos",
) -> pd.DataFrame:
    """
    Use mygene.info to map gene identifiers to another type.

    Parameters
    ----------
    query: iterable
        a list or list-like of gene identifiers
    tax_id: str or int
        Target genome taxonomy id
    field : str, optional
        Target identifier to map the query genes to. Valid fields
        are: ensembl.gene, entrezgene, symbol, name, refseq, entrezgene. Note that
        refseq will return the protein refseq_id by default, use `refseq.translation.rna` to
        return the RNA refseq_id. Currently, mapping to Ensembl transcript ids is
        not supported.

    Returns
    -------
    pandas.DataFrame
        mapped gene annotation.
    """
    field, _ = _parse_mygene_input(field)

    logger.info("Querying mygene.info...")
    query = list(query)
    query_len = len(query)
    batch_size = 1000  # same as mygene.info internal batch size
    it = range(0, query_len, batch_size)
    if query_len > batch_size:
        it = tqdm(it, unit=f"{batch_size} queries")

    result = pd.DataFrame()
    mg = mygene.MyGeneInfo()
    for i in it:
        _result = mg.querymany(
            query[i : i + batch_size],  # noqa
            scopes="symbol,name,ensembl.gene,entrezgene,ensembl.transcript,"
            "ensembl,accession.protein,accession.rna,other_names,alias",
            field=field,
            species=tax_id,
            as_dataframe=True,
            verbose=False,
        )
        result = pd.concat((result, _result))

    if "notfound" in result and result.shape[1] == 1:
        logger.warning("No matching genes found")

    return result


def _parse_mygene_input(field, product=None):
    if product is not None:
        product = product.lower()
        if product not in ["rna", "protein"]:
            raise ValueError("Argument 'product' should be either 'rna' or 'protein'.")

    field = field.lower()
    # see mg.MyGeneInfo().get_fields()
    allowed_fields = [
        "ensembl.gene",
        "entrezgene",
        "symbol",
        "name",
        "refseq",
        "genomic_pos",
    ]
    if product is None:
        allowed_fields.append("refseq.translation.rna")
    if field not in allowed_fields:
        raise ValueError(
            f"Argument 'field' should be either in {', '.join(allowed_fields)}"
        )

    return field, product


def _filter_query(query: pd.DataFrame) -> pd.DataFrame:
    """per queried gene, keep the best matching, non-NaN, mapping"""
    if "notfound" in query:
        query = query[query.notfound.isnull()]  # drop unmatched genes
        query = query.drop(columns=["notfound"])

    query = query.dropna()

    if "query" in query:
        query = query.groupby("query").first()  # already sorted by mapping score
        query = query.drop(columns=["_id", "_score"])
    return query
