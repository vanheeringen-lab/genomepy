import re
from typing import List

from loguru import logger

from genomepy.exceptions import GenomeDownloadError
from genomepy.files import get_file_info
from genomepy.online import read_url
from genomepy.providers.base import BaseProvider
from genomepy.utils import get_localname


class UrlProvider(BaseProvider):
    """
    URL genome provider.

    Simply download a genome directly through a url.
    """

    name = "URL"
    _cli_install_options = {
        "to_annotation": {
            "long": "to-annotation",
            "help": "link to the annotation file, required if this is not in the same directory as the fasta file",
            "default": None,
        },
    }

    def __init__(self):
        self.genomes = {}

    @staticmethod
    def ping():
        return True

    def genome_taxid(self, name):
        return

    def assembly_accession(self, name):
        return

    def search(self, term):
        """return an empty generator,
        same as if no genomes were found at the other providers"""
        yield from ()

    def _genome_info_tuple(self, name):
        return tuple()

    def _check_name(self, name):
        """check if genome name can be found for provider"""
        return name

    def get_genome_download_link(self, url, mask=None, **kwargs):
        return url

    def get_annotation_download_link(self, name: str, **kwargs) -> str:
        """
        Return a functional annotation download link.

        Parameters
        ----------
        name : str
            genome name
        **kwargs: dict, optional:
            to_annotation : direct URL to the gene annotation

        Returns
        -------
        str
            http/ftp link

        Raises
        ------
        GenomeDownloadError
            if no functional link was found
        """
        link = kwargs.get("to_annotation")
        if link:
            ext = get_file_info(link)[0]
            if ext not in [".gtf", ".gff", ".gff3", ".bed"]:
                raise TypeError(
                    "Only (gzipped) gtf, gff and bed files are supported.\n"
                )
            return link

        links = self.get_annotation_download_links(name)
        if links:
            return links[0]

        raise GenomeDownloadError(
            f"No gene annotations found for {get_localname(name)}.\n"
        )

    def get_annotation_download_links(self, name: str, **kwargs) -> List[str]:
        """
        Retrieve functioning gene annotation download link(s).

        If provided, check if the annotation url links to a supported file type (gtf/gff3/bed).
        Else try to find an annotation in the same location as the genome url.

        Parameters
        ----------
        name : str
            genome name

        Returns
        -------
        list
            http/ftp link(s)
        """
        # name = url to genome
        return search_url_for_annotations(name)


def search_url_for_annotations(url: str) -> list:
    """Attempts to find gtf or gff3 files in the same location as the genome url"""
    name = get_localname(url)

    urldir = url[: url.rfind("/")]
    logger.info(
        "You have requested the gene annotation to be downloaded. "
        "Genomepy will check the remote directory: "
        f"{urldir} for annotation files..."
    )

    # try to find a GTF or GFF3 file
    dirty_list = read_url(urldir).split("\n")
    fnames = fuzzy_annotation_search(name, dirty_list)
    links = [urldir + "/" + fname for fname in fnames]
    if not links:
        logger.warning(
            "Could not parse the remote directory. "
            "Please supply a URL using --URL-to-annotation.\n"
        )
    return links


def fuzzy_annotation_search(search_name, search_list):
    """Returns all files containing both name and an annotation extension"""
    hits = []
    for ext in ["gtf", "gff"]:
        # .*? = non greedy filler. 3? = optional 3 (for gff3). (\.gz)? = optional .gz
        expr = f"{search_name}.*?\.{ext}3?(\.gz)?"  # noqa: W605
        for line in search_list:
            hit = re.search(expr, line, flags=re.IGNORECASE)
            if hit:
                hits.append(hit[0])
    return hits
