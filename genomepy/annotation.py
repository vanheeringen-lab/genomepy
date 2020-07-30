import sys
from typing import Iterable, Optional, Union, TextIO

from loguru import logger
import pandas as pd

from genomepy.provider import ProviderBase
from genomepy import Genome


##@cached
def load_mapping(
    to: str, provider: Optional[str] = None, fmt: Optional[str] = "dataframe"
) -> Union[pd.DataFrame, dict]:
    """Load chromosome mapping from one version of assembly to another

    Parameters
    ----------
    to: str
        Either a local genomepy genome name, a remote genome name or an assembly accession (GCA*)

    provider: str, optional
        If an accession is specified the provider needs to be specified (UCSC, Ensembl or NCBI)

    fmt: str, optional
        Format for the chromosome mapping, either "dataframe" or "dict".

    Returns
    -------
    pandas.DataFrame or dict
        Chromosome mapping.
    """
    if fmt.lower() not in ["dataframe", "dict"]:
        raise ValueError("Invalid format, should be 'dataframe' or 'dict'")

    logger.info("Loading chromosome mapping.")
    if to.startswith("GCA"):
        if provider is None:
            raise ValueError("Need a provider: NCBI, UCSC or Ensembl")
        asm_acc = to
    else:
        try:
            genome = Genome(to)
            logger.info("Using local genome information")
            asm_acc = genome.assembly_accession
            if provider is None:
                provider = genome.provider
            asm_report = genome.assembly_report
        except Exception:
            logger.info("Searching remote genome information")
            result = [row for row in ProviderBase.search(to, provider=provider)]
            if len(result) > 1:
                p = [row[1] for row in result]
                raise ValueError(
                    f"More than one result, need one of these providers: {', '.join(p)}"
                )
            if provider is None:
                provider = result[0][1]
            asm_acc = result[0][2]
            asm_report = ProviderBase.create(provider).download_assembly_report(asm_acc)

    logger.info(f"Assembly {asm_acc}, provider {provider}")

    if provider not in ["UCSC", "NCBI", "Ensembl"]:
        logger.error(f"Can't map to provider {provider}")
        return None

    asm_report.loc[
        asm_report["Sequence-Role"] != "assembled-molecule", "Assigned-Molecule"
    ] = "na"

    asm_report["Ensembl-Name"] = asm_report["Sequence-Name"]
    asm_report.loc[
        asm_report["Sequence-Role"] != "assembled-molecule", "Ensembl-Name"
    ] = asm_report.loc[
        asm_report["Sequence-Role"] != "assembled-molecule", "GenBank-Accn"
    ]

    if provider == "NCBI":
        logger.info("Mapping to NCBI sequence names")
        id_column = "Sequence-Name"
    elif provider == "UCSC":
        logger.info("Mapping to UCSC sequence names")
        id_column = "UCSC-style-name"
        ucsc_ids = asm_report[id_column].unique()
        if len(ucsc_ids) == 1 and ucsc_ids[0] == "na":
            raise ValueError("UCSC style names not available for this assembly")
    elif provider == "Ensembl":
        logger.info("Mapping to Ensembl sequence names")
        id_column = "Ensembl-Name"

    mapping = asm_report[
        [
            "Sequence-Name",
            "UCSC-style-name",
            "Assigned-Molecule",
            "GenBank-Accn",
            "Ensembl-Name",
        ]
    ]
    mapping = pd.melt(mapping, id_vars=[id_column])
    mapping = mapping[mapping["value"] != "na"]
    mapping = mapping.drop_duplicates().set_index("value")[[id_column]]
    mapping.columns = ["chrom"]
    if fmt == "dataframe":
        return mapping
    if fmt == "dict":
        return mapping["chrom"].to_dict()


def map_bedfile(
    infile: str,
    to: str,
    path_or_buf: Optional[Union[None, TextIO, str]] = None,
    provider: Optional[str] = None,
    keepgoing: Optional[bool] = False,
):
    """Map BED file to specific assembly version.

    Parameters
    ----------
    infile: str
        Filename of BED file.

    to: str
        Either a local genomepy genome name, a remote genome name or an assembly accession (GCA*)

    path_or_buf : str or filehandle, optional
        Filehandle or name of output file. If None (default), then sys.stdout will be used.

    provider: str, optional
        If an accession is specified the provider needs to be specified (UCSC, Ensembl or NCBI)

    keepgoing: bool, optional
        Ignore errors if chromosome mappings are not found.
    """
    if path_or_buf is None:
        fout = sys.stdout
    elif isinstance(path_or_buf, str):
        fout = open(path_or_buf, "w")
    else:
        fout = path_or_buf

    map_dict = load_mapping(to, provider, fmt="dict")
    seen = {}
    for line in open(infile):
        vals = line.strip().split("\t")
        chrom = vals[0]
        new_chrom = map_dict.get(chrom, None)
        if new_chrom is None:
            if keepgoing:
                if not seen.get(chrom):
                    logger.warning(
                        f"Skipping all features for {chrom}, not found in mapping."
                    )
                    seen[chrom] = 1
                continue
            else:
                raise ValueError(f"No mapping found for {chrom}")
        print("\t".join([new_chrom] + vals[1:]), file=fout)

    if isinstance(path_or_buf, str):
        fout.close()


# @cached
def local_genome_coords(genes: Iterable[str], genome: Genome) -> pd.DataFrame:
    """Retrieve gene location from local annotation.

    Parameters
    ----------
    genes : Iterable
        List of gene names or gene identifiers such as ensembl_id.
    genome : str
        Genome name

    Returns
    -------
    pandas.DataFrame with gene annotation.
    """
    gene_list = list(genes)
    gene_info = genome.gene_annotation(genome)
    if gene_info is None:
        return

    gene_info = gene_info[["chrom", "start", "end", "name", "strand"]].set_index("name")
    gene_info = gene_info.loc[gene_list]

    # If we find more than half of the genes we assume this worked.
    if gene_info.shape[0] >= 0.5 * len(gene_list):
        return gene_info.reset_index()[["chrom", "start", "end", "name", "strand"]]


# @cached
def genome_coords(genes: Iterable[str], genome: str) -> pd.DataFrame:
    """Retrieve genomic coordinates of a set of genes.

    If the annotation is not present locally, then mygene.info is used.
    All genome assemblies that are present in the latest version of
    Ensembl are supported by mygene.info.

    Parameters
    ----------
    genes : Iterable
        List of gene names or gene identifiers such as ensembl_id.
    genome : str
        Genome name

    Returns
    -------
    pandas.DataFrame with gene annotation.
    """
    g = Genome(genome)

    # First try to find the genes in the local annotation installed by genomepy.
    gene_info = local_genome_coords(genes, genome)
    if gene_info is not None:
        return gene_info

    # Genes are not identified locally.
    # Retrieve the gene information using the mygene.info API
    logger.info(f"No local matching genes found for {g.name}, trying mygene.info")
    result = g._query_mygene(genes)

    if g.provider == "Ensembl":
        result = result.rename(columns={"genomic_pos.chr": "chrom"})
    else:
        # Ensembl, UCSC and NCBI chromosome names can all be different :-/
        logger.info("Local genome is not an Ensembl genome.")
        mapping = load_mapping(g.name)
        result = result.join(mapping, on="genomic_pos.chr")
        result = result.dropna(subset=["chrom"])

    # Convert genomic positions from string to integer
    result["genomic_pos.start"] = result["genomic_pos.start"].astype(int)
    result["genomic_pos.end"] = result["genomic_pos.end"].astype(int)

    # For each gene use the match with the highest score
    result = result.reset_index().sort_values("_score").groupby("query").last()

    # Map the Ensembl 1/-1 strand to +/- strand
    strand_df = pd.DataFrame({"ens_strand": [-1, 1], "strand": ["-", "+"]}).set_index(
        "ens_strand"
    )
    result = result.join(strand_df, on="genomic_pos.strand")

    # Select the correct columns and name them
    result = result.reset_index()[
        ["chrom", "genomic_pos.start", "genomic_pos.end", "query", "strand"]
    ]
    result.columns = [["chrom", "start", "end", "name", "strand"]]

    return result
