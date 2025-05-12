import re

import requests
from loguru import logger

from genomepy.caching import cache_exp_genomes, cache_exp_other, disk_cache, lock
from genomepy.exceptions import GenomeDownloadError
from genomepy.online import check_url, retry
from genomepy.providers.base import BaseProvider
from genomepy.utils import safe


class EnsemblProvider(BaseProvider):
    """
    Ensembl genome provider.

    Will search both ensembl.org and ensemblgenomes.org.
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
        api_online = bool(check_url("https://rest.ensembl.org/info/ping?"))
        vertebrate_url_online = bool(check_url("http://ftp.ensembl.org"))
        other_url_online = bool(check_url("http://ftp.ensemblgenomes.org"))
        return api_online and vertebrate_url_online and other_url_online

    def _genome_info_tuple(self, name, size=False):
        """tuple with assembly metadata"""
        accession = self.assembly_accession(name)
        taxid = self.genome_taxid(name)
        annotations = bool(self.annotation_links(name))
        species = self.genomes[name].get("scientific_name")
        other = self.genomes[name].get("genebuild")
        if size:
            length = self.genomes[name]["base_count"]
            return name, accession, taxid, annotations, species, length, other
        return name, accession, taxid, annotations, species, other

    def get_version(self, name: str, version=None) -> int:
        """
        Retrieve the latest Ensembl or EnsemblGenomes release version,
        or check if the requested release version exists.
        """
        division, is_vertebrate = self.get_division(name)
        if version is None:
            latest_version = self.get_release(is_vertebrate)
            return latest_version

        if not str(version).isdecimal():
            raise TypeError("Version must be a number")
        version = int(version)

        all_versions = self.get_releases(is_vertebrate)
        ensembl = f"Ensembl{'' if is_vertebrate else 'Genomes'}"
        if version not in all_versions:
            raise ValueError(
                f"{ensembl} release version {version} "
                f"not found. Available versions: {all_versions}"
            )

        releases = self.releases_with_assembly(name)
        if version not in releases:
            raise FileNotFoundError(
                f"{name} not found on {ensembl} release {version}. "
                f"Available on release versions: {releases}"
            )
        return version

    def get_division(self, name: str):
        """Retrieve the division of a genome."""
        genome = self.genomes[safe(name)]
        division = str(genome["division"]).lower().replace("ensembl", "")
        if division == "bacteria":
            raise NotImplementedError("Bacteria from Ensembl not supported.")

        is_vertebrate = division == "vertebrates"
        return division, is_vertebrate

    @disk_cache.memoize(
        expire=cache_exp_other, tag="get_release-ensembl", ignore={"self"}
    )
    def get_release(self, is_vertebrate: bool) -> int:
        """Retrieve current Ensembl or EnsemblGenomes release version."""
        ext = "/info/data/?" if is_vertebrate else "/info/eg_version?"
        ret = retry(request_json, 3, self._url, ext)
        return int(ret["releases"][0] if is_vertebrate else ret["version"])

    @disk_cache.memoize(expire=cache_exp_other, ignore={"self"})
    def get_releases(self, is_vertebrate: bool):
        """Retrieve all Ensembl or EnsemblGenomes release versions."""
        url = "http://ftp.ensemblgenomes.org/pub?"
        if is_vertebrate:
            url = "http://ftp.ensembl.org/pub?"
        ret = retry(requests.get, 3, url)
        # sort releases new to old
        releases = sorted(
            [int(i) for i in re.findall(r'"release-(\d+)/"', ret.text)],
            reverse=True,
        )
        if is_vertebrate:
            # ignore immature releases
            releases = [r for r in releases if r > 46]

        # exclude pre-releases (returns 403: forbidden)
        latest_release = self.get_release(is_vertebrate)
        releases = [r for r in releases if r <= latest_release]
        return releases

    @lock
    @disk_cache.memoize(
        expire=cache_exp_other, tag="get_releases-ensembl", ignore={"self"}
    )
    def releases_with_assembly(self, name: str):
        """List all Ensembl or EnsemblGenomes release versions with the specified genome."""
        genome = self.genomes[safe(name)]
        lwr_name = genome["name"]
        asm_name = re.sub(r"\.p\d+$", "", safe(genome["assembly_name"]))
        division, is_vertebrate = self.get_division(name)

        # all releases with the genome fasta
        releases = self.get_releases(is_vertebrate)
        releases_with_assembly = []
        for release in releases:
            url = f"http://ftp.ensemblgenomes.org/pub/release-{release}/{division}/fasta/{lwr_name}/dna/"
            if is_vertebrate:
                url = f"https://ftp.ensembl.org/pub/release-{release}/fasta/{lwr_name}/dna/"
            ret = retry(requests.get, 3, url)
            if asm_name in ret.text:  # 404 error has text too, so this always works
                releases_with_assembly.append(release)
            else:
                break
        return releases_with_assembly

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
        division, is_vertebrate = self.get_division(name)

        # base directory of the genome
        ftp = "http://ftp.ensemblgenomes.org"
        if is_vertebrate:
            ftp = "http://ftp.ensembl.org"
        version = self.get_version(name, kwargs.get("version"))
        div_path = "" if is_vertebrate else f"/{division}"
        lwr_name = genome["name"]
        ftp_directory = f"{ftp}/pub/release-{version}{div_path}/fasta/{lwr_name}/dna"

        # this assembly has its own directory
        if name == "GRCh37":
            ftp_directory = genome["genome"].format(version)

        # specific fasta file
        cap_name = lwr_name.capitalize()
        asm_name = re.sub(r"\.p\d+$", "", safe(genome["assembly_name"]))
        mask_lvl = {"soft": "_sm", "hard": "_rm", "none": ""}[mask]
        asm_lvl = "toplevel" if kwargs.get("toplevel") else "primary_assembly"
        version_tag = "" if version > 30 else f".{version}"

        ftp_file = f"{cap_name}.{asm_name}{version_tag}.dna{mask_lvl}.{asm_lvl}.fa.gz"

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
            version : Ensembl version to use. By default, the latest version is used

        Returns
        -------
        list
            http link(s)
        """
        genome = self.genomes[safe(name)]
        division, is_vertebrate = self.get_division(name)

        # base directory of the genome
        ftp = "http://ftp.ensemblgenomes.org"
        if is_vertebrate:
            ftp = "http://ftp.ensembl.org"
        version = self.get_version(name, kwargs.get("version"))
        div_path = "" if is_vertebrate else f"/{division}"
        lwr_name = genome["name"]
        ftp_directory = f"{ftp}/pub/release-{version}{div_path}/gtf/{lwr_name}"

        # specific gtf file
        cap_name = lwr_name.capitalize()
        asm_name = re.sub(r"\.p\d+$", "", safe(genome["assembly_name"]))

        ftp_file = f"{cap_name}.{asm_name}.{version}.gtf.gz"

        # combine
        link = f"{ftp_directory}/{ftp_file}"
        if name == "GRCh37":
            link = genome["annotation"].format(version)
        return [link] if check_url(link, max_tries=2) else []


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
        "name": "Homo_sapiens",
        "scientific_name": "Homo sapiens",
        "url_name": "human",
        "assembly_name": "GRCh37",
        "division": "vertebrates",
        "base_count": "3137144693",
        "display_name": "",
        "genebuild": "",
    }
    return genomes


@lock
@disk_cache.memoize(expire=cache_exp_genomes, tag="get_genomes-ensembl")
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

        # filter summaries to these keys (to reduce the size of the cached data)
        summary_keys_to_keep = [
            "assembly_name",
            "assembly_accession",
            "taxonomy_id",
            "name",
            "scientific_name",
            "url_name",
            "display_name",
            "genebuild",
            "division",
            "base_count",
        ]
        for genome in division_genomes:
            if "_gca_" in genome["name"]:
                continue  # ~1600 mislabeled protists and fungi
            name = safe(genome["assembly_name"])
            genomes[name] = {k: genome[k] for k in summary_keys_to_keep}

    genomes = add_grch37(genomes)
    return genomes
