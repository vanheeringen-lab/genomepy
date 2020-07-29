import os
import sys
from typing import Tuple, Iterable, Optional, Union

from loguru import logger
import mygene
from genomepy.provider import ProviderBase #, cached
from genomepy import Genome, search
import pandas as pd
import numpy as np


logger.remove()
logger.add(
    sys.stderr,
    format="<green>{time:YYYY-MM-DD at HH:mm:ss}</green> <bold>|</bold> <blue>{level}</blue> <bold>|</bold> {message}",
    level="INFO",
)


# @cached
def ensembl_genome_info(genome_name: str) -> Tuple[str, str, str]:
    """Return Ensembl genome information for a local genome managed by genomepy.

    Parameters
    ----------
    genome_name : str
        Name of local genome.

    Returns
    -------
    (str, str, str)
        Ensembl name, accession, taxonomy_id
    """
    # Fast lookup for some common queries
    common_names = {
        "danRer11": "GRCz11",
        "hg38": "GRCh38.p13",
        "mm10": "GRCm38.p6",
        "dm6": "BDGP6.28",
    }
    if genome_name in common_names:
        search_term = common_names[genome_name]
    else:
        try:
            genome = Genome(genome_name)
            search_term = genome.tax_id
        except FileNotFoundError:
            logger.info(f"Genome {genome_name} not installed locally")
            p = ProviderBase.create("Ensembl")
            for name, *_rest in p.search(genome_name):
                if name == genome_name:
                    logger.info(
                        f"It can be downloaded from Ensembl: genomepy install {name} Ensembl --annotation"
                    )
                    return None
            return None

    # search Ensembl by taxonomy_id or by specific Ensembl name (if we know it)
    p = ProviderBase.create("Ensembl")
    logger.info(f"searching for {search_term}")
    results = list(p.search(str(search_term)))
    if len(results) == 0:
        logger.warn(f"Could not find a genome for this species in Ensembl.")
    
    for name, accession, species, tax_id, *rest in results:
        # Check if the assembly_id of the current Ensembl genome is the same as the
        # local genome. If it is identical, we can correctly assume that the genomes
        # sequences are identical.
        # For the genomes in the lookup table, we already know they match.
        if common_names[genome_name] == name or accession == genome.assembly_accession:
            return name, accession, tax_id

    logger.warn(f"Could not find a matching assembly version in the current release of Ensembl.")

    return None


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
            asm_report = genome.assembly_reportr
        except Exception:
            logger.info("Searching remote genome information")
            result = [row for row in search(to, provider=provider)]
            if len(result) > 1:
                p = [row[1].decode() for row in result]
                raise ValueError(
                    f"More than one result, need one of these providers: {', '.join(p)}"
                )
            if provider is None:
                provider = result[0][1].decode()
            asm_acc = result[0][2].decode()
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

    provider: str, optional
        If an accession is specified the provider needs to be specified (UCSC, Ensembl or NCBI)

    keepgoing: bool, optional
        Ignore errors if chromosome mappings are not found.
    """
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
        print("\t".join([new_chrom] + vals[1:]))


# @cached
def map_gene_dataframe(
    df: pd.DataFrame, genome: str, gene_field: str, product: str = "protein"
) -> pd.DataFrame:
    """Use mygene.info to map identifiers

    If the identifier can't be mapped, it will be dropped from the resulting
    annotation.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with gene information. Should contain at least a "name" column.

    genome : str
        Genome name for mygene.info.

    gene_field : str
        Identifier for gene annotation. Uses mygene.info to map ids. Valid fields
        are: ensembl.gene, entrezgene, symbol, name, refseq, entrezgene. Note that
        refseq will return the protein refseq_id by default, use `product="rna"` to
        return the RNA refseq_id. Currently, mapping to Ensembl transcript ids is
        not supported.

    product : str, optional
        Either "protein" or "rna". Only used when `gene_field="refseq"`

    Returns
    -------
    pandas.DataFrame with mapped gene annotation.
    """
    cols = df.columns
    df["split_id"] = df["name"].str.split(r"\.", expand=True)[0]
    genes = df["split_id"].tolist()
    result = _query_mygene(genome, genes, fields=gene_field)
    if result is None:
        logger.error("Could not map using mygene.info")
    df = df.join(result, on="split_id")

    # Only in case of RefSeq annotation the product needs to be specified.
    if gene_field == "refseq":
        gene_field = f"{gene_field}.translation.{product}"

    # Get rid of extra columns from query
    df["name"] = df[gene_field]
    df = df[cols].dropna()

    return df


# @cached
def gene_annotation(
    genome: str, gene_field: Optional[str] = None, product: str = "protein"
) -> pd.DataFrame:
    """Retrieve gene location from local annotation.

    You can use mygene.info to map identifiers on the fly by specifying
    `gene_field`. If the identifier can't be mapped, it will be dropped
    from the resulting annotation.

    Returns a DataFrame with gene annotation in bed12 format.

    Parameters
    ----------
    genome : str
        Genome name

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
    product = product.lower()
    gene_field = gene_field.lower()

    if product not in ["rna", "protein"]:
        raise ValueError(f"Argument product should be either 'rna' or 'protein'")

    bed12_fields = [
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

    g = Genome(genome)

    if gene_field is not None:
        if gene_field == "refseq":
            anno_file = f"{genome}.{gene_field}.{product}.annotation.bed"
        else:
            anno_file = f"{genome}.{gene_field}.annotation.bed"

        target_anno = os.path.join(os.path.dirname(g.filename), anno_file)
        if os.path.exists(target_anno):
            return pd.read_csv(target_anno, sep="\t", names=bed12_fields, dtype={"chrom":"string", "start":np.uint32, "end":np.uint32})

    for anno_file in [f"{genome}.annotation.bed.gz", f"{genome}.annotation.bed"]:
        bed = os.path.join(os.path.dirname(g.filename), anno_file)
        if os.path.exists(bed):
            break
        else:
            bed = None

    if bed is None:
        logger.info(f"No annotation file found for genome {genome}")
        return

    df = pd.read_csv(bed, sep="\t", names=bed12_fields)

    # Optionally use mygene.info to map gene/transcript ids
    if gene_field is not None:
        logger.info("Mapping gene identifiers using mygene.info")
        logger.info("This can take a long time, but results will be saved for faster access next time!")
        df = map_gene_dataframe(df, genome, gene_field=gene_field, product=product)
        
        # Save mapping results for quicker access next time
        df.to_csv(target_anno, sep="\t", index=False, header=False)

    return df


# @cached
def local_genome_coords(genes: Iterable[str], genome: str) -> pd.DataFrame:
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
    gene_info = gene_annotation(genome)
    if gene_info is None:
        return

    gene_info = gene_info[["chrom", "start", "end", "name", "strand"]].set_index("name")
    gene_info = gene_info.loc[gene_list]

    # If we find more than half of the genes we assume this worked.
    if gene_info.shape[0] >= 0.5 * len(gene_list):
        return gene_info.reset_index()[["chrom", "start", "end", "name", "strand"]]


def _query_mygene(genome: str, query: Iterable[str], fields: str = "genomic_pos"):
    # mygene.info only queries the most recent version of the Ensembl database
    # We can only safely continue if the local genome matched the Ensembl genome.
    # Even if the local genome was installed via Ensembl, we still need to check
    # if it is the same version
    result = ensembl_genome_info(genome)
    if result is None:
        return None

    # Run the actual query
    g = Genome(genome)
    logger.info("Querying mygene.info...")
    mg = mygene.MyGeneInfo()
    result = mg.querymany(
        query,
        scopes="symbol,name,ensembl.gene,entrezgene,ensembl.transcript,ensembl",
        fields=fields,
        species=g.tax_id,
        as_dataframe=True,
        verbose=False,
    )
    logger.info("Done")
    if "notfound" in result and result.shape[1] == 1:
        logger.error("No matching genes found")
        sys.exit()

    return result


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
    # First try to find the genes in the local annotation installed by genomepy.
    gene_info = local_genome_coords(genes, genome)
    if gene_info is not None:
        return gene_info

    # Genes are not identified locally.
    # Retrieve the gene information using the mygene.info API
    logger.info(f"No local matching genes found for {genome}, trying mygene.info")
    result = _query_mygene(genome, genes)

    g = Genome(genome)
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
