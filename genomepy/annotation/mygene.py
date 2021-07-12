from typing import Iterable, Optional, Tuple, Union

import mygene
import pandas as pd
from loguru import logger
from tqdm import tqdm

from genomepy.annotation.utils import _check_property, _parse_annot
from genomepy.caching import memory
from genomepy.files import read_readme
from genomepy.providers import nearest_assembly
from genomepy.utils import safe


def map_genes(
    self,
    gene_field: str,
    product: str = "protein",
    annot: Union[str, pd.DataFrame] = "bed",
) -> Optional[pd.DataFrame]:
    """
    Uses mygene.info to map gene identifiers on the fly by specifying
    `gene_field`. If the identifier can't be mapped, it will be dropped
    from the resulting annotation.

    Returns the dataframe with remapped "name" column.

    Parameters
    ----------
    self: Annotation class instance

    annot: str or pd.Dataframe
        Annotation dataframe to map (a pandas dataframe or "bed").
        Is mapped to a column named "name" (required).

    gene_field : str, optional
        Identifier for gene annotation. Uses mygene.info to map ids. Valid fields
        are: ensembl.gene, entrezgene, symbol, name, refseq, entrezgene. Note that
        refseq will return the protein refseq_id by default, use `product="rna"` to
        return the RNA refseq_id. Currently, mapping to Ensembl transcript ids is
        not supported.

    product : str, optional
        Either "protein" or "rna". Only used when `gene_field="refseq"`

    Returns
    -------
    pandas.DataFrame with gene annotation.
    """
    to, product = parse_mygene_input(gene_field, product)
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

    result = _query_mygene(self, genes, fields=to)
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
    field : str, optional
        Target identifier to map the query genes to. Valid fields
        are: ensembl.gene, entrezgene, symbol, name, refseq, entrezgene. Note that
        refseq will return the protein refseq_id by default, use `refseq.translation.rna` to
        return the RNA refseq_id. Currently, mapping to Ensembl transcript ids is
        not supported.
    batch_size: int, optional
        Controls batch size for the REST API.

    Returns
    -------
    pandas.DataFrame with mapped gene annotation.
    """
    field, _ = parse_mygene_input(field)

    logger.info("Querying mygene.info...")
    query = list(set(query))
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
            fields=field,
            species=tax_id,
            as_dataframe=True,
            verbose=False,
        )
        result = pd.concat((result, _result))

    if "notfound" in result and result.shape[1] == 1:
        logger.warning("No matching genes found")

    return result


def parse_mygene_input(gene_field, product=None):
    if product:
        product = product.lower()
        if product not in ["rna", "protein"]:
            raise ValueError("Argument product should be either 'rna' or 'protein'.")

    gene_field = gene_field.lower()
    allowed_fields = [
        "ensembl.gene",
        "entrezgene",
        "symbol",
        "name",
        "refseq",
        "entrezgene",
    ]
    if product is None:
        allowed_fields.append("refseq.translation.rna")
    if gene_field not in allowed_fields:
        raise ValueError(
            "Argument product should be either 'ensembl.gene', "
            "'entrezgene,' 'symbol', 'name', 'refseq' or 'entrezgene'."
        )

    return gene_field, product


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
    search_result = nearest_assembly(asm_acc, "Ensembl")
    if search_result is None:
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
    ensembl_info = ensembl_genome_info(self)
    if ensembl_info is None:
        return pd.DataFrame()

    tax_id = ensembl_info[2]
    if not str(tax_id).isdigit():
        raise ValueError("No taxomoy ID found")

    return query_mygene(query, tax_id, fields, batch_size)


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
