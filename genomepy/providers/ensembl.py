import re

import requests
from loguru import logger

from genomepy.caching import cache, goldfish_cache
from genomepy.exceptions import GenomeDownloadError
from genomepy.online import check_url, retry
from genomepy.providers.base import BaseProvider
from genomepy.utils import safe


class EnsemblProvider(BaseProvider):
    """
    Ensembl genome provider.

    Will search both ensembl.org as well as ensemblgenomes.org.
    The bacteria division is not yet supported.
    """

    name = "Ensembl"
    accession_fields = ["assembly_accession"]
    taxid_fields = ["taxonomy_id"]
    description_fields = [
        "name",
        "scientific_name",
        "url_name",
        "display_name",
    ]
    _cli_install_options = {
        "toplevel": {
            "long": "toplevel",
            "help": "always download toplevel-genome",
            "flag_value": True,
        },
        "version": {
            "long": "version",
            "help": "select release version",
            "type": int,
            "default": None,
        },
    }
    _url = "https://rest.ensembl.org/"

    def __init__(self):
        self._provider_status()
        # Populate on init, so that methods can be cached
        self.genomes = get_genomes(self._url)

    @staticmethod
    def ping():
        """Can the provider be reached?"""
        return bool(check_url("https://rest.ensembl.org/info/ping?"))

    def _genome_info_tuple(self, name):
        """tuple with assembly metadata"""
        accession = self.assembly_accession(name)
        taxid = self.genome_taxid(name)
        annotations = bool(self.annotation_links(name))
        species = self.genomes[name].get("scientific_name")
        other = self.genomes[name].get("genebuild")
        return name, accession, taxid, annotations, species, other

    @goldfish_cache(ignore=["self"])
    def get_version(self, vertebrates=False, set_version=None):
        """Retrieve current version from Ensembl FTP."""
        if set_version:
            return str(set_version)
        ext = "/info/data/?" if vertebrates else "/info/eg_version?"
        ret = retry(request_json, 3, self._url, ext)
        releases = ret["releases"] if vertebrates else [ret["version"]]
        return str(max(releases))

    def get_genome_download_link(self, name, mask="soft", **kwargs):
        """
        Return http link to the genome sequence

        Parameters
        ----------
        name : str
            Genome name. Current implementation will fail if exact
            name is not found.

        mask : str , optional
            Masking level. Options: soft, hard or none. Default is soft.

        Returns
        ------
        str with the http download link.
        """
        genome = self.genomes[safe(name)]
        division, is_vertebrate = get_division(genome)

        # base directory of the genome
        ftp = "http://ftp.ensemblgenomes.org"
        if is_vertebrate:
            ftp = "http://ftp.ensembl.org"
        version = self.get_version(is_vertebrate, kwargs.get("version"))
        div_path = "" if is_vertebrate else f"/{division}"
        lwr_name = genome["url_name"].lower()

        ftp_directory = f"{ftp}/pub/release-{version}{div_path}/fasta/{lwr_name}/dna"
        # this assembly has its own directory
        if name == "GRCh37":
            ftp_directory = genome["genome"].format(version)

        # specific fasta file
        cap_name = genome["url_name"].capitalize()
        asm_name = re.sub(r"\.p\d+$", "", safe(genome["assembly_name"]))
        mask_lvl = {"soft": "_sm", "hard": "_rm", "none": ""}[mask]
        asm_lvl = "toplevel" if kwargs.get("toplevel") else "primary_assembly"

        ftp_file = f"{cap_name}.{asm_name}.dna{mask_lvl}.{asm_lvl}.fa.gz"

        # combine
        link = f"{ftp_directory}/{ftp_file}"
        if check_url(link, 2):
            return link

        # primary assemblies do not always exist
        if asm_lvl == "primary_assembly":
            link = link.replace("primary_assembly", "toplevel")
            if check_url(link, 2):
                return link

        raise GenomeDownloadError(
            f"Could not download genome {name} from {self.name}.\n"
            "URL is broken. Select another genome or provider.\n"
            f"Broken URL: {link}"
        )

    def get_annotation_download_links(self, name, **kwargs):
        """
        Retrieve functioning gene annotation download link(s).

        Parameters
        ----------
        name : str
            genome name
        **kwargs: dict, optional:
            version : Ensembl version to use. By default the latest version is used

        Returns
        -------
        list
            http link(s)
        """
        genome = self.genomes[safe(name)]
        division, is_vertebrate = get_division(genome)

        # base directory of the genome
        ftp = "http://ftp.ensemblgenomes.org"
        if is_vertebrate:
            ftp = "http://ftp.ensembl.org"
        version = self.get_version(is_vertebrate, kwargs.get("version"))
        div_path = "" if is_vertebrate else f"/{division}"
        lwr_name = genome["url_name"].lower()

        ftp_directory = f"{ftp}/pub/release-{version}{div_path}/gtf/{lwr_name}"

        # specific gtf file
        cap_name = genome["url_name"].capitalize()
        asm_name = re.sub(r"\.p\d+$", "", safe(genome["assembly_name"]))

        ftp_file = f"{cap_name}.{asm_name}.{version}.gtf.gz"

        # combine
        link = f"{ftp_directory}/{ftp_file}"
        if name == "GRCh37":
            link = genome["annotation"].format(version)
        return [link] if check_url(link, max_tries=2) else []


def get_division(genome: dict):
    """Retrieve the division of a genome."""
    division = genome["division"].lower().replace("ensembl", "")
    if division == "bacteria":
        raise NotImplementedError("Bacteria from Ensembl not supported.")

    is_vertebrate = bool(division == "vertebrates")
    return division, is_vertebrate


def request_json(rest_url, ext):
    """Make a REST request and return as json."""
    if rest_url.endswith("/") and ext.startswith("/"):
        ext = ext[1:]

    r = requests.get(rest_url + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()

    return r.json()


def add_grch37(genomes):
    genomes["GRCh37"] = {
        "genome": (
            "http://ftp.ensembl.org/pub/grch37/release-{}/fasta/homo_sapiens/dna"
        ),
        "annotation": (
            "http://ftp.ensembl.org/pub/grch37/release-{}/gtf/"
            "homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"
        ),
        "assembly_accession": "GCA_000001405.14",
        "taxonomy_id": 9606,
        "name": "human",
        "scientific_name": "Homo sapiens",
        "url_name": "Homo_sapiens",
        "assembly_name": "GRCh37",
        "division": "vertebrates",
        "display_name": "",
        "genebuild": "",
    }
    return genomes


@cache
def get_genomes(rest_url):
    logger.info("Downloading assembly summaries from Ensembl")

    genomes = {}
    divisions = retry(request_json, 3, rest_url, "info/divisions?")
    for division in divisions:
        if division == "EnsemblBacteria":
            continue
        division_genomes = retry(
            request_json, 3, rest_url, f"info/genomes/division/{division}?"
        )
        for genome in division_genomes:
            genomes[safe(genome["assembly_name"])] = genome

    genomes = add_grch37(genomes)
    return genomes
