import os
import re
import urllib.error
from typing import List
from urllib.request import urlopen

import pandas as pd
from loguru import logger
from tqdm.auto import tqdm

from genomepy.caching import cache
from genomepy.exceptions import GenomeDownloadError
from genomepy.online import check_url, read_url
from genomepy.providers.base import BaseProvider
from genomepy.utils import mkdir_p, safe

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


class NcbiProvider(BaseProvider):
    """
    NCBI genome provider.

    Uses the assembly reports page to search and list genomes.
    """

    name = "NCBI"
    accession_fields = ["assembly_accession", "gbrs_paired_asm"]
    taxid_fields = ["species_taxid", "taxid"]
    description_fields = [
        "submitter",
        "organism_name",
        "assembly_accession",
        "gbrs_paired_asm",
        "paired_asm_comp",
    ]
    _cli_install_options = {}
    _url = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/"

    def __init__(self):
        self._provider_status()
        # Populate on init, so that methods can be cached
        self.genomes = get_genomes(self._url)

    @staticmethod
    def ping():
        """Can the provider be reached?"""
        return bool(check_url("https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/"))

    def _genome_info_tuple(self, name):
        """tuple with assembly metadata"""
        accession = self.assembly_accession(name)
        taxid = self.genome_taxid(name)
        annotations = bool(self.annotation_links(name))
        species = self.genomes[name].get("organism_name")
        other = self.genomes[name].get("submitter")
        return name, accession, taxid, annotations, species, other

    def get_genome_download_link(self, name, mask="soft", **kwargs):
        """
        Return NCBI ftp link to top-level genome sequence

        Parameters
        ----------
        name : str
            Genome name. Current implementation will fail if exact
            name is not found.

        mask : str , optional
            Masking level. Options: soft, hard or none. Default is soft.

        Returns
        ------
        str with the http/ftp download link.
        """
        # only soft masked genomes available. can be (un)masked in _post_process_download
        link = self._ftp_or_html_link(name, file_suffix="_genomic.fna.gz")

        if link:
            return link

        raise GenomeDownloadError(
            f"Could not download genome {name} from {self.name}.\n"
            "URL is broken. Select another genome or provider.\n"
            f"Broken URL: {link}"
        )

    def _post_process_download(self, name, fname, out_dir, mask="soft"):
        """
        Replace accessions with sequence names in fasta file.

        Applies masking.

        Parameters
        ----------
        name : str
            NCBI genome name

        fname : str
            file path to the genome fasta

        out_dir : str
            output directory

        mask : str , optional
            masking level: soft/hard/none, default=soft
        """
        # (down)load the assembly report
        asm_fname = os.path.join(out_dir, "assembly_report.txt")
        if not os.path.exists(asm_fname):
            # save assembly report
            download_assembly_report(self.assembly_accession(name), asm_fname)
        asm_report = pd.read_csv(
            asm_fname, sep="\t", comment="#", header=None, dtype=str
        )

        # create mapping of chromosome accessions to chromosome names
        tr = asm_report.set_index(6)[0].to_dict()

        # mask sequence if required
        if mask == "soft":

            def mask_cmd(txt):
                return txt

        elif mask == "hard":
            logger.info("NCBI genomes are softmasked by default. Hard masking...")
            lower_case_nucleotides = re.compile("[actg]")

            def mask_cmd(txt):
                return lower_case_nucleotides.sub("N", txt)

        else:
            logger.info("NCBI genomes are softmasked by default. Unmasking...")

            def mask_cmd(txt):
                return txt.upper()

        # apply mapping and masking
        old_fname = os.path.join(
            os.path.dirname(fname), f"original_{os.path.basename(fname)}"
        )
        os.rename(fname, old_fname)
        with open(old_fname) as old, open(fname, "w") as new:
            for line in tqdm(
                old, desc="Processing NCBI Fasta", unit_scale=1, unit=" lines"
            ):
                if line[0] == ">":
                    desc = line[1:].strip()
                    name = desc.split(" ")[0]
                    new.write(f">{tr.get(name, name)} {desc}\n")
                else:
                    new.write(mask_cmd(line))

    def get_annotation_download_links(self, name: str, **kwargs) -> List[str]:
        """
        Retrieve functioning gene annotation download link(s).

        Parse and test the link to the NCBI annotation file.

        Parameters
        ----------
        name : str
            genome name

        Returns
        -------
        list
            http/ftp link(s)
        """
        link = self._ftp_or_html_link(name, file_suffix="_genomic.gff.gz")
        return [link] if link else []

    def _ftp_or_html_link(self, name, file_suffix, skip_check=False):
        """
        NCBI's files are accessible over FTP and HTTPS
        Try HTTPS first and return the first functioning link
        """
        genome = self.genomes[safe(name)]
        ftp_link = genome["ftp_path"]
        html_link = ftp_link.replace("ftp://", "https://")
        for link in [html_link, ftp_link]:
            link += "/" + link.split("/")[-1] + file_suffix

            if skip_check or check_url(link, max_tries=2, timeout=10):
                return link


@cache
def get_genomes(assembly_url):
    """Parse genomes from assembly summary txt files."""
    logger.info("Downloading assembly summaries from NCBI, this will take a while...")

    def load_summary(url):
        """
        lazy loading of the url so we can parse while downloading
        """
        for row in urlopen(url):
            yield row

    genomes = {}
    # order is important as asm_name can repeat (overwriting the older name)
    names = [
        "assembly_summary_genbank_historical.txt",
        "assembly_summary_refseq_historical.txt",
        "assembly_summary_genbank.txt",
        "assembly_summary_refseq.txt",
    ]
    for fname in names:
        lines = load_summary(f"{assembly_url}/{fname}")
        _ = next(lines)  # line 0 = comment
        header = (
            next(lines).decode("utf-8").strip("# ").strip("\n").split("\t")
        )  # line 1 = header
        for line in tqdm(lines, desc=fname[17:-4], unit_scale=1, unit=" genomes"):
            line = line.decode("utf-8").strip("\n").split("\t")
            if line[19] != "na":  # ftp_path must exist
                name = safe(line[15])  # overwrites older asm_names
                genomes[name] = dict(zip(header, line))
    return genomes


def download_assembly_report(acc: str, fname: str = None):
    """
    Retrieve the NCBI assembly report.

    Returns the assembly_report as a pandas DataFrame if fname is not specified.

    Parameters
    ----------
    acc : str
        Assembly accession (GCA or GCF)
    fname : str, optional
        Save assembly_report to this filename.

    Returns
    -------
    pandas.DataFrame
        NCBI assembly report.
    """
    msg = "Could not download the assembly report from NCBI. "
    if not isinstance(acc, str) or not acc.startswith(("GCA", "GCF")):
        logger.warning(msg)
        return None
    assembly_report = _assembly_report_url(acc)
    if assembly_report is None:
        logger.warning(msg + f"Assembly accession '{acc}' not found.")
        return None
    asm_report = pd.read_csv(
        assembly_report, sep="\t", comment="#", names=ASM_FORMAT, dtype=str
    )

    if fname:
        mkdir_p(os.path.dirname(fname))
        asm_report.to_csv(fname, sep="\t", index=False)
    else:
        return asm_report


def _assembly_report_url(acc: str) -> str or None:
    ftp_basedir = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
    asm_dir = f"{ftp_basedir}/{acc[0:3]}/{acc[4:7]}/{acc[7:10]}/{acc[10:13]}/"
    try:
        text = read_url(asm_dir)
    except urllib.error.HTTPError:
        return None

    # try to find exact accession patch level
    hits = re.findall(acc + '_.*/"', text)  # len 0-1
    if len(hits) == 0:
        # patch level not found, pick the closest instead
        hits = re.findall(acc + '.*/"', text)
        available_accessions = ["_".join(h.split("_")[0:2]) for h in hits]
        closest = _closest_patch_lvl(acc, available_accessions)
        hits = [h for h in hits if closest + "_" in h]  # len 0-1
    if len(hits) == 0:
        return None

    name = hits[0][:-2]
    assembly_report = asm_dir + name + "/" + name + "_assembly_report.txt"
    return assembly_report


def _closest_patch_lvl(reference: str, targets: list) -> str:
    ref_patch = _patch_lvl(reference)
    closest = ""
    closest_dist = 999
    for target in targets:
        tgt_patch = _patch_lvl(target)
        dist = abs(tgt_patch - 0.1 - ref_patch)  # tiebreaker: newer > older patches
        if dist < closest_dist:
            closest = target
            closest_dist = dist
    return closest


def _patch_lvl(acc: str) -> int:
    """str(GCA_000000000.6) -> int(6)"""
    return int(acc.split(".")[1]) if "." in acc else 0
