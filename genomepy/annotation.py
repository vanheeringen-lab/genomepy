import os
import sys
from typing import Tuple, Iterable, Optional

from loguru import logger
import mygene
from genomepy.provider import ProviderBase, cached
from genomepy import Genome
import pandas as pd


logger.remove()
logger.add(
    sys.stderr,
    format="<green>{time:YYYY-MM-DD at HH:mm:ss}</green> <bold>|</bold> <blue>{level}</blue> <bold>|</bold> {message}",
    level="INFO",
)


@cached
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
        "hg38": "GRCh38",
        "mm10": "GRCm38",
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
    name, accession, species, tax_id, *rest = [row for row in p.search(search_term)][0]

    # Check if the assembly_id of the current Ensembl genome is the same as the
    # local genome. If it is identical, we can correctly assume that the genomes
    # sequences are identical.
    # For the genomes in the lookup table, we already know they match.
    if genome_name in common_names or accession == genome.assembly_accession:
        return name, accession, tax_id
    else:
        print(f"Could not find a matching genome in Ensembl")
        return None


@cached
def ncbi_assembly_report(asm_acc: str) -> pd.DataFrame:
    """Retrieve the NCBI assembly report as a DataFrame.

    Parameters
    ----------
    asm_acc : str
        Assembly accession (GCA or GCF)

    Returns
    -------
    pandas.DataFrame
        NCBI assembly report.
    """
    p = ProviderBase.create("NCBI")
    ncbi_search = list(p.search(asm_acc))
    if len(ncbi_search) > 1:
        raise Exception("More than one genome for accession")
    else:
        ncbi_name = ncbi_search[0][0].replace(" ", "_")

    # NCBI FTP location of assembly report
    logger.info(f"Found NCBI assembly {asm_acc} with name {ncbi_name}")
    assembly_report = (
        f"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{asm_acc[0:3]}/"
        + f"{asm_acc[4:7]}/{asm_acc[7:10]}/{asm_acc[10:13]}/"
        + f"{asm_acc}_{ncbi_name}/{asm_acc}_{ncbi_name}_assembly_report.txt"
    )

    logger.info(f"Downloading {assembly_report}")
    header = [
        "Sequence-Name",
        "Sequence-Role",
        "Assigned-Molecule",
        "Assigned-Molecule-Location/Type",
        "GenBank-Accn",
        "Relationship",
        "RefSeq-Accn",
        "Assembly-Unit",
        "Sequence-Length",
        "UCSC-style-name",
    ]
    asm_report = pd.read_csv(assembly_report, sep="\t", comment="#", names=header)
    return asm_report


@cached
def load_mapping(genome_name):
    logger.info("Loading chromosome mapping.")
    genome = Genome(genome_name)
    asm_acc = genome.assembly_accession

    if genome.provider not in ["UCSC", "NCBI"]:
        logger.error(f"Can't map to provider {genome.provider}")
        return None

    asm_report = ncbi_assembly_report(asm_acc)
    asm_report.loc[
        asm_report["Sequence-Role"] != "assembled-molecule", "Assigned-Molecule"
    ] = "na"

    mapping = asm_report[
        ["Sequence-Name", "UCSC-style-name", "Assigned-Molecule", "GenBank-Accn"]
    ]

    if genome.provider == "NCBI":
        logger.info("Mapping to NCBI sequence names")
        id_column = "Sequence-Name"
    elif genome.provider == "UCSC":
        logger.info("Mapping to UCSC sequence names")
        id_column = "UCSC-style-name"
    mapping = pd.melt(mapping, id_vars=[id_column])
    mapping = mapping[mapping["value"] != "na"]
    mapping = mapping.drop_duplicates().set_index("value")[[id_column]]
    mapping.columns = ["chrom"]
    return mapping


@cached
def map_gene_dataframe(df: pd.DataFrame, genome: str, gene_field: str, product: str = "protein") -> pd.DataFrame:
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
    df = df.join(result, on="split_id")

    # Only in case of RefSeq annotation the product needs to be specified.
    if gene_field == "refseq":
        gene_field = f"{gene_field}.translation.{product}"

    # Get rid of extra columns from query
    df["name"] = df[gene_field]
    df = df[cols].dropna()

    return df


@cached
def gene_annotation(
    genome: str, gene_field: Optional[str] = None, product: str = "protein"
) -> pd.DataFrame:
    """Retrieve gene location from local annotation.

    You can use mygene.info to map identifiers on the fly by specifying
    `gene_field`. If the identifier can't be mapped, it will be dropped
    from the resulting annotation.

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
    if product not in ["rna", "protein"]:
        raise ValueError(f"Argument product should be either 'rna' or 'protein'")

    g = Genome(genome)
    for anno_file in [f"{genome}.annotation.bed.gz", f"{genome}.annotation.bed"]:
        bed = os.path.join(os.path.dirname(g.filename), anno_file)
        if os.path.exists(bed):
            break
        else:
            bed = None

    if bed is None:
        logger.info(f"No annotation file found for genome {genome}")
        return

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
    df = pd.read_table(bed, names=bed12_fields)

    # Optionally use mygene.info to map gene/transcript ids
    if gene_field is not None:
        df = map_gene_dataframe(df, genome, gene_field=gene_field, product=product)

    return df


@cached
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


@cached
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