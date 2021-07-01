import os
import re
from typing import Iterator

import requests
from loguru import logger

from genomepy.caching import cache
from genomepy.exceptions import GenomeDownloadError
from genomepy.online import check_url, read_url
from genomepy.providers.base import BaseProvider
from genomepy.providers.ncbi import NcbiProvider
from genomepy.utils import lower


class UcscProvider(BaseProvider):
    """
    UCSC genome provider.

    The UCSC API REST server is used to search and list genomes.
    """

    name = "UCSC"
    accession_fields = []
    taxid_fields = ["taxId"]
    description_fields = ["description", "scientificName"]
    _cli_install_options = {
        "ucsc_annotation_type": {
            "long": "annotation",
            "help": "specify annotation to download: UCSC, Ensembl, NCBI_refseq or UCSC_refseq",
            "default": None,
        },
    }
    _url = "http://hgdownload.soe.ucsc.edu/goldenPath"

    def __init__(self):
        self._provider_status()
        # Populate on init, so that methods can be cached
        self.genomes = get_genomes("http://api.genome.ucsc.edu/list/ucscGenomes")

    @staticmethod
    def ping():
        """Can the provider be reached?"""
        return bool(check_url("http://hgdownload.soe.ucsc.edu/goldenPath"))

    def _search_accession(self, term: str) -> Iterator[str]:
        """
        UCSC does not store assembly accessions.
        This function searches NCBI (most genomes + stable accession IDs),
        then uses the NCBI accession search results for a UCSC text search.

        Parameters
        ----------
        term : str
            Assembly accession, GCA_/GCF_....

        Yields
        ------
        genome names
        """
        # NCBI provides a consistent assembly accession. This can be used to
        # retrieve the species, and then search for that.
        p = NcbiProvider()
        ncbi_genomes = list(p._search_accession(term))

        # remove superstrings (keep GRCh38, not GRCh38.p1 to GRCh38.p13)
        unique_ncbi_genomes = []
        for i in ncbi_genomes:
            if sum([j in i for j in ncbi_genomes]) == 1:
                unique_ncbi_genomes.append(i)

        # add NCBI organism names to search terms
        organism_names = [
            p.genomes[name]["organism_name"] for name in unique_ncbi_genomes
        ]
        terms = list(set(unique_ncbi_genomes + organism_names))

        # search with NCBI results in the given provider
        for name, metadata in self.genomes.items():
            for term in terms:
                term = lower(term)
                if term in lower(name) or any(
                    [term in lower(metadata[f]) for f in self.description_fields]
                ):
                    yield name
                    break  # max one hit per genome

    @cache(ignore=["self"])
    def assembly_accession(self, name: str) -> str:
        """Return the assembly accession (GCA_/GCF_....) for a genome.

        UCSC does not serve the assembly accession through the REST API.
        Therefore, the readme.html is scanned for an assembly accession. If it is
        not found, the linked NCBI assembly page will be checked.

        Parameters
        ----------
        name: str
            genome name

        Returns
        ------
        str
            Assembly accession.
        """
        try:
            ucsc_url = (
                "https://hgdownload.soe.ucsc.edu/" + self.genomes[name]["htmlPath"]
            )
            text = read_url(ucsc_url)
        except UnicodeDecodeError:
            return "na"

        # example accessions: GCA_000004335.1 (ailMel1)
        # regex: GC[AF]_ = GCA_ or GCF_, \d = digit, \. = period
        accession_regex = re.compile(r"GC[AF]_\d{9}\.\d+")
        match = accession_regex.search(text)
        if match:
            return match.group(0)

        # Search for an assembly link at NCBI
        match = re.search(r"https?://www.ncbi.nlm.nih.gov/assembly/\d+", text)
        if match:
            ncbi_url = match.group(0)
            text = read_url(ncbi_url)

            # retrieve valid assembly accessions.
            # contains additional info, such as '(latest)' or '(suppressed)'. Unused for now.
            valid_accessions = re.findall(r"assembly accession:.*?GC[AF]_.*?<", text)
            text = " ".join(valid_accessions)
            match = accession_regex.search(text)
            if match:
                return match.group(0)

    def _annot_types(self, name):
        links = self.annotation_links(name)
        if links:
            annot_types = ["knownGene", "ensGene", "ncbiRefSeq", "refGene"]
            annotations_found = [False, False, False, False]
            for n, annot in enumerate(annot_types):
                for link in links:
                    if annot in link:
                        annotations_found[n] = True
            return annotations_found

    def _genome_info_tuple(self, name):
        """tuple with assembly metadata"""
        accession = self.assembly_accession(name)
        taxid = self.genome_taxid(name)
        annotations = self._annot_types(name)
        species = self.genomes[name].get("scientificName")
        other = self.genomes[name].get("description")
        return name, accession, taxid, annotations, species, other

    def get_genome_download_link(self, name, mask="soft", **kwargs):
        """
        Return UCSC http link to genome sequence

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
        ucsc_url = self._url + "/{0}/bigZips/chromFa.tar.gz"
        ucsc_url_masked = self._url + "/{0}/bigZips/chromFaMasked.tar.gz"
        alt_ucsc_url = self._url + "/{0}/bigZips/{0}.fa.gz"
        alt_ucsc_url_masked = self._url + "/{0}/bigZips/{0}.fa.masked.gz"

        # soft masked genomes. can be unmasked in _post _process_download
        urls = [ucsc_url, alt_ucsc_url]
        if mask == "hard":
            urls = [ucsc_url_masked, alt_ucsc_url_masked]

        for genome_url in urls:
            link = genome_url.format(name)

            if check_url(link, 2):
                return link

        raise GenomeDownloadError(
            f"Could not download genome {name} from {self.name}.\n"
            "URLs are broken. Select another genome or provider.\n"
            f"Broken URLs: {', '.join([url.format(name) for url in urls])}"
        )

    @staticmethod
    def _post_process_download(name, localname, out_dir, mask="soft"):  # noqa
        """
        Unmask a softmasked genome if required

        Parameters
        ----------
        name : str
            unused for the UCSC function

        localname : str
            Custom name for your genome

        out_dir : str
            Output directory

        mask : str , optional
            masking level: soft/hard/none, default=soft
        """
        if mask != "none":
            return

        logger.info("UCSC genomes are softmasked by default. Unmasking...")

        fa = os.path.join(out_dir, f"{localname}.fa")
        old_fa = os.path.join(out_dir, f"old_{localname}.fa")
        os.rename(fa, old_fa)
        with open(old_fa) as old, open(fa, "w") as new:
            for line in old:
                if line.startswith(">"):
                    new.write(line)
                else:
                    new.write(line.upper())

    def get_annotation_download_links(self, name, **kwargs):
        """
        Parse and test the link to the UCSC annotation file.

        Will check UCSC, Ensembl, NCBI RefSeq and UCSC RefSeq annotation, respectively.
        More info on the annotation file on: https://genome.ucsc.edu/FAQ/FAQgenes.html#whatdo

        Parameters
        ----------
        name : str
            Genome name
        """
        gtf_url = f"http://hgdownload.soe.ucsc.edu/goldenPath/{name}/bigZips/genes/"
        txt_url = f"http://hgdownload.cse.ucsc.edu/goldenPath/{name}/database/"

        # download gtf format if possible, txt format if not
        gtfs_exists = check_url(gtf_url, 2)
        base_url = gtf_url + name + "." if gtfs_exists else txt_url
        base_ext = ".gtf.gz" if gtfs_exists else ".txt.gz"

        links = []
        for annot in ["knownGene", "ensGene", "ncbiRefSeq", "refGene"]:
            link = base_url + annot + base_ext
            if check_url(link, max_tries=2):
                links.append(link)
        if links:
            return links

    def get_annotation_download_link(self, name: str, **kwargs) -> str:
        """select a functional annotation download link from a list of links"""
        links = self.annotation_links(name)
        if links:
            # user specified annotation type
            annot_type = kwargs.get("ucsc_annotation_type", "").lower()
            if annot_type:
                annot_files = {
                    "ucsc": "knownGene",
                    "ensembl": "ensGene",
                    "ncbi_refseq": "ncbiRefSeq",
                    "ucsc_refseq": "refGene",
                }
                if annot_type not in annot_files:
                    raise ValueError(
                        f"UCSC annotation type '{annot_type}' not recognized."
                    )

                annot = annot_files[annot_type]
                links = [link for link in links if annot in link]
                if links:
                    return links[0]
                raise FileNotFoundError(
                    f"Assembly {name} does not possess the {annot_type} '{annot}' gene annotation."
                )

            # first working link
            return links[0]
        raise FileNotFoundError(f"No gene annotations found for {name} on {self.name}")


@cache
def get_genomes(rest_url):
    logger.info("Downloading assembly summaries from UCSC")

    r = requests.get(rest_url, headers={"Content-Type": "application/json"})
    if not r.ok:
        r.raise_for_status()
    ucsc_json = r.json()
    genomes = ucsc_json["ucscGenomes"]
    return genomes
