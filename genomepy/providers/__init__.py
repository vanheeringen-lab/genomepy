"""Provider class, modules & related functions"""
import os
from typing import List, Optional

import pandas as pd
from loguru import logger

from genomepy.files import read_readme
from genomepy.providers.base import ASM_FORMAT
from genomepy.providers.ensembl import EnsemblProvider
from genomepy.providers.gencode import GencodeProvider
from genomepy.providers.local import LocalProvider
from genomepy.providers.ncbi import NcbiProvider
from genomepy.providers.ucsc import UcscProvider
from genomepy.providers.url import UrlProvider
from genomepy.utils import get_genomes_dir, safe

__all__ = [
    "Provider",
    "list_online_providers",
    "list_providers",
    "search",
    "map_locations",
    "download_assembly_report",
    "nearest_assembly",
    "online_providers",
]

PROVIDERS = {
    "gencode": GencodeProvider,
    "ensembl": EnsemblProvider,
    "ucsc": UcscProvider,
    "ncbi": NcbiProvider,
    "local": LocalProvider,
    "url": UrlProvider,
}


def create(name: str):
    """
    Create a provider based on the provider name.

    Adds additional method(s) to the instance.

    Parameters
    ----------
    name : str
        Name of the provider (e.g. UCSC, Ensembl, ...)

    Returns
    -------
    provider
        Provider instance
    """
    name = name.lower()
    if name not in PROVIDERS:
        options = "', '".join(list_providers())
        raise ValueError(f"Unknown provider '{name}'. Options: '{options}'")

    p = PROVIDERS[name]
    p.download_assembly_report = staticmethod(download_assembly_report)
    return p()


def list_providers():
    """
    List of providers genomepy supports

    Returns
    -------
    list
        names of providers
    """
    return [p.name for p in PROVIDERS.values()]


def list_online_providers():
    """
    List of providers genomepy supports *that are online right now*.

    Returns
    -------
    list
        names of online providers
    """
    return [p.name for p in PROVIDERS.values() if p.ping()]


def online_providers(provider: str = None):
    """
    Check if the provider can be reached, or any provider if none is specified.

    Parameters
    ----------
    provider : str, optional
        Only try to yield the specified provider.

    Yields
    ------
    provider
        Provider instances
    """
    for provider in [provider] if provider else list_providers():
        try:
            yield create(provider)
        except ConnectionError as e:
            logger.warning(str(e))


def search(term, provider: str = None):
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
        Only search the specified provider (faster).

    Yields
    ------
    list
        genome name, provider and metadata
    """
    term = safe(str(term))
    for p in online_providers(provider):
        for row in p.search(term):
            ret = list(row[:1]) + [p.name] + list(row[1:])
            yield ret


def _closest_patch_lvl(reference, targets):
    ref_patch = int(reference.split(".")[1]) if "." in reference else 0
    nearest = [999, []]
    for target in targets:
        tgt_patch = int(target.split(".")[1]) if "." in target else 0
        distance = abs(tgt_patch + 0.1 - ref_patch)  # tiebreaker: newer > older patches
        if distance == nearest[0]:
            nearest[1].append(target)
        if distance < nearest[0]:
            nearest = [distance, [target]]
    return nearest[1]


def _best_accession(reference: str, targets: list):
    """Return the nearest accession ID from a list of IDs"""
    if len(targets) == 1:
        return targets[0]

    # GCA/GCF
    matching_prefix = [t for t in targets if t.startswith(reference[0:3])]
    if len(matching_prefix) > 0:
        targets = matching_prefix

    # patch levels
    # e.g. GCA_000002035.4 & GCA_000002035.3
    targets = _closest_patch_lvl(reference, targets)

    if len(targets) > 1:
        logger.info(
            f"Multiple matching accession numbers found. Selecting the closest ({targets[0]})."
        )
    return targets[0]


def _best_search_result(asm_acc: str, results: List[list]) -> list:
    """Return the best search result based on accession IDs."""
    results = [res for res in results if res[2] is not None]

    if len(results) > 1:
        accessions = [res[2] for res in results]
        bes_acc = _best_accession(asm_acc, accessions)
        results = [res for res in results if res[2] == bes_acc]

    if len(results) > 0:
        return results[0]


def nearest_assembly(asm_acc: str, provider: str) -> list:
    """
    Return the search result of the assembly nearest to
    the given accession ID, in the specified provider.
    """
    if not asm_acc.startswith(("GCA_", "GCF_")):
        raise ValueError("asm_acc must be an assembly accession ID.")

    search_results = list(search(asm_acc, provider=provider))
    best_result = _best_search_result(asm_acc, search_results)
    if best_result is None:
        logger.warning(f"No assembly similar to {asm_acc} on {provider}.")
    return best_result


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
    search_result = nearest_assembly(asm_acc, "NCBI")
    if search_result is None:
        return
    ncbi_acc = search_result[2]
    ncbi_name = search_result[0]

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
        logger.warning("Cannot map without an assembly report.")
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
        logger.warning("UCSC style names not available for this assembly.")
        return

    mapping = asm_report[[f"{frm_provider}_name", f"{to_provider}_name"]]
    mapping = mapping.dropna().drop_duplicates().set_index(f"{frm_provider}_name")
    return mapping


class Provider:
    """Find & instantiate providers"""

    list = staticmethod(list_providers)
    online = staticmethod(list_online_providers)
    create = staticmethod(create)
