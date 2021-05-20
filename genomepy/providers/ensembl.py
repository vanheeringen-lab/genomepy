import re

from loguru import logger
import requests

from genomepy.exceptions import GenomeDownloadError
from genomepy.online import check_url, retry
from genomepy.providers.base import BaseProvider, cache, goldfish_cache
from genomepy.utils import safe


class EnsemblProvider(BaseProvider):
    """
    Ensembl genome provider.

    Will search both ensembl.org as well as ensemblgenomes.org.
    The bacteria division is not yet supported.
    """

    accession_fields = ["assembly_accession"]
    taxid_fields = ["taxonomy_id"]
    description_fields = [
        "name",
        "scientific_name",
        "url_name",
        "display_name",
    ]
    provider_specific_install_options = {
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

    def __init__(self):
        self.name = "Ensembl"
        self.rest_url = "https://rest.ensembl.org/"
        self.provider_status(self.rest_url + "info/ping?", max_tries=2)
        # Populate on init, so that methods can be cached
        self.genomes = self._get_genomes(self.rest_url)

    @cache(ignore=["self"])
    def _get_genomes(self, rest_url):
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
        return genomes

    def _genome_info_tuple(self, name):
        """tuple with assembly metadata"""
        accession = self.assembly_accession(name)
        taxid = self.genome_taxid(name)
        annotations = bool(self.annotation_links(name))
        species = self.genomes[name].get("scientific_name")
        other = self.genomes[name].get("genebuild")
        return name, accession, taxid, annotations, species, other

    @goldfish_cache(ignore=["self", "rest_url"])
    def get_version(self, rest_url, vertebrates=False):
        """Retrieve current version from Ensembl FTP."""
        ext = "/info/data/?" if vertebrates else "/info/eg_version?"
        ret = retry(request_json, 3, rest_url, ext)
        releases = ret["releases"] if vertebrates else [ret["version"]]
        return str(max(releases))

    def get_genome_download_link(self, name, mask="soft", **kwargs):
        """
        Return Ensembl http or ftp link to the genome sequence

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
        genome = self.genomes[safe(name)]

        # parse the division
        division = genome["division"].lower().replace("ensembl", "")
        if division == "bacteria":
            raise NotImplementedError("bacteria from ensembl not yet supported")

        ftp_site = "ftp://ftp.ensemblgenomes.org/pub"
        if division == "vertebrates":
            ftp_site = "ftp://ftp.ensembl.org/pub"

        # Ensembl release version
        version = kwargs.get("version")
        if version is None:
            version = self.get_version(self.rest_url, division == "vertebrates")

        # division dependent url format
        ftp_dir = "{}/release-{}/fasta/{}/dna".format(
            division, version, genome["url_name"].lower()
        )
        if division == "vertebrates":
            ftp_dir = "release-{}/fasta/{}/dna".format(
                version, genome["url_name"].lower()
            )
        url = f"{ftp_site}/{ftp_dir}"

        # masking and assembly level
        def get_url(level="toplevel"):
            masks = {"soft": "dna_sm.{}", "hard": "dna_rm.{}", "none": "dna.{}"}
            pattern = masks[mask].format(level)

            asm_url = "{}/{}.{}.{}.fa.gz".format(
                url,
                genome["url_name"].capitalize(),
                re.sub(r"\.p\d+$", "", safe(genome["assembly_name"])),
                pattern,
            )
            return asm_url

        # try to get the (much smaller) primary assembly,
        # unless specified otherwise
        link = get_url("primary_assembly")
        if kwargs.get("toplevel") or not check_url(link, 2):
            link = get_url()

        if check_url(link, 2):
            return link

        raise GenomeDownloadError(
            f"Could not download genome {name} from {self.name}.\n"
            "URL is broken. Select another genome or provider.\n"
            f"Broken URL: {link}"
        )

    def get_annotation_download_links(self, name, **kwargs):
        """
        Parse and test the link to the Ensembl annotation file.

        Parameters
        ----------
        name : str
            Genome name
        kwargs: dict , optional:
            Provider specific options.

            version : int , optional
                Ensembl version. By default the latest version is used.
        """
        genome = self.genomes[safe(name)]
        division = genome["division"].lower().replace("ensembl", "")

        ftp_site = "ftp://ftp.ensemblgenomes.org/pub"
        if division == "vertebrates":
            ftp_site = "ftp://ftp.ensembl.org/pub"

        # Ensembl release version
        version = kwargs.get("version")
        if version is None:
            version = self.get_version(self.rest_url, division == "vertebrates")

        if division != "vertebrates":
            ftp_site += f"/{division}"

        # Get the GTF URL
        base_url = ftp_site + "/release-{}/gtf/{}/{}.{}.{}.gtf.gz"
        safe_name = re.sub(r"\.p\d+$", "", name)
        link = base_url.format(
            version,
            genome["url_name"].lower(),
            genome["url_name"].capitalize(),
            safe_name,
            version,
        )

        if check_url(link, max_tries=2):
            return [link]


def request_json(rest_url, ext):
    """Make a REST request and return as json."""
    if rest_url.endswith("/") and ext.startswith("/"):
        ext = ext[1:]

    r = requests.get(rest_url + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()

    return r.json()
