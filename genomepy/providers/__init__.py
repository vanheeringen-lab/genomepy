import os
from typing import Optional

import pandas as pd
from loguru import logger

from genomepy.files import read_readme
from genomepy.providers.ensembl import EnsemblProvider
from genomepy.providers.ncbi import NcbiProvider
from genomepy.providers.ucsc import UcscProvider
from genomepy.providers.url import UrlProvider
from genomepy.utils import best_search_result, get_genomes_dir, safe

ASM_FORMAT = [
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

PROVIDERS = {
    "ensembl": EnsemblProvider,
    "ucsc": UcscProvider,
    "ncbi": NcbiProvider,
    "url": UrlProvider,
}


def create(name: str):
    """
    Create a provider based on the provider name.

    Parameters
    ----------
    name : str
        Name of the provider (eg. UCSC, Ensembl, ...)

    Returns
    -------
    provider :
        Provider instance.
    """
    name = name.lower()
    if name not in PROVIDERS:
        logger.error(f"Unknown provider. Options: {', '.join(list_providers())}")
        return

    p = PROVIDERS[name]()
    return p


def list_providers():
    """List available providers."""
    return list(PROVIDERS.keys())


def online_providers(provider: str = None):
    """
    Check if the provider can be reached, or any provider if none is specified.
    Return online provider(s) as objects.
    """
    for provider in [provider] if provider else list_providers():
        try:
            yield create(provider)
        except ConnectionError as e:
            logger.warning(str(e))


def search_all(term, provider: str = None, encode: bool = False):
    """
    Search for a genome.

    If provider is specified, search only that specific provider, else
    search all providers. Both the name and description are used for the
    search. Search term is case-insensitive.

    Parameters
    ----------
    term : str
        Search term, case-insensitive.
    provider : str , optional
        Provider name
    encode : bool, optional
        Encode return strings.

    Yields
    ------
    tuple
        genome information (name/identifier and description)
    """
    term = safe(str(term))
    for p in online_providers(provider):
        for row in p.search(term):
            ret = list(row[:1]) + [p.name] + list(row[1:])
            if encode:
                ret = [x.encode("utf-8") for x in ret]
            yield ret


def nearest_ncbi_assembly(asm_acc):
    """
    return the NCBI accession ID and NCBI assembly name nearest to
    the given accession ID.
    """
    ncbi_search = list(NcbiProvider().search(asm_acc))
    search_result = best_search_result(asm_acc, ncbi_search, acc_idx=1)
    if len(search_result) == 0:
        logger.warning(f"Could not download an NCBI assembly report for {asm_acc}")
        return None, None
    ncbi_acc = search_result[1]
    ncbi_name = search_result[0]
    return ncbi_acc, ncbi_name


def download_assembly_report(asm_acc: str, fname: str = None):
    """
    Retrieve the NCBI assembly report.

    Returns the assembly_report as a pandas DataFrame if fname is not specified.

    Parameters
    ----------
    asm_acc : str
        Assembly accession (GCA or GCF)
    fname : str, optional
        Save assembly_report to this filename.

    Returns
    -------
    pandas.DataFrame
        NCBI assembly report.
    """
    ncbi_acc, ncbi_name = nearest_ncbi_assembly(asm_acc)
    if ncbi_acc is None:
        return

    # NCBI FTP location of assembly report
    assembly_report = (
        f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{ncbi_acc[0:3]}/"
        + f"{ncbi_acc[4:7]}/{ncbi_acc[7:10]}/{ncbi_acc[10:13]}/"
        + f"{ncbi_acc}_{ncbi_name}/{ncbi_acc}_{ncbi_name}_assembly_report.txt"
    )
    asm_report = pd.read_csv(assembly_report, sep="\t", comment="#", names=ASM_FORMAT)

    if fname:
        asm_report.to_csv(fname, sep="\t", index=False)
    else:
        return asm_report


def map_locations(
    frm: str, to: str, genomes_dir: Optional[str] = None
) -> Optional[pd.DataFrame]:
    """
    Load chromosome mapping from one assembly to another using the
    NCBI assembly reports.

    Parameters
    ----------
    frm: str
        A local genomepy genome name
    to: str
        target provider (UCSC, Ensembl or NCBI)
    genomes_dir: str, optional
        The genomes directory to look for the genomes.
        Will search the default genomes_dir if left blank.

    Returns
    -------
    pandas.DataFrame
        Chromosome mapping.
    """
    to_provider = to.lower()
    if to_provider not in ["ucsc", "ncbi", "ensembl"]:
        raise ValueError(f"Genomepy can only map to NCBI, UCSC or Ensembl, not '{to}'.")

    genomes_dir = get_genomes_dir(genomes_dir)
    frm_readme = os.path.join(genomes_dir, frm, "README.txt")
    frm_asm_report = os.path.join(genomes_dir, frm, "assembly_report.txt")
    if not os.path.exists(frm_readme):
        raise FileNotFoundError(f"Cannot find {frm} in {genomes_dir}.")

    metadata, _ = read_readme(frm_readme)
    frm_provider = metadata.get("provider").lower()
    if frm_provider == to_provider:
        logger.warning(f"You are attempting to map {frm} from {to} to {to}.")
        return

    asm_acc = metadata.get("assembly_accession")
    if not os.path.exists(frm_asm_report):
        download_assembly_report(asm_acc, frm_asm_report)
    if not os.path.exists(frm_asm_report):
        logger.error("Cannot map without an assembly report")
        return

    asm_report = pd.read_csv(frm_asm_report, sep="\t", comment="#")
    asm_report["ensembl_name"] = asm_report["Sequence-Name"]
    asm_report["ncbi_name"] = asm_report["Sequence-Name"]
    asm_report["ucsc_name"] = asm_report["UCSC-style-name"]

    # for Ensembl, use GenBank names for the scaffolds
    asm_report.loc[
        asm_report["Sequence-Role"] != "assembled-molecule", "ensembl_name"
    ] = asm_report.loc[
        asm_report["Sequence-Role"] != "assembled-molecule", "GenBank-Accn"
    ]

    if "ucsc" in [frm_provider, to_provider] and list(
        asm_report["ucsc_name"].unique()
    ) == ["na"]:
        logger.error("UCSC style names not available for this assembly")
        return

    mapping = asm_report[[f"{frm_provider}_name", f"{to_provider}_name"]]
    mapping = mapping.dropna().drop_duplicates().set_index(f"{frm_provider}_name")
    return mapping
