import os
import re

from genomepy.exceptions import GenomeDownloadError
from genomepy.files import get_file_info
from genomepy.providers.base import BaseProvider
from genomepy.utils import cleanpath, get_genomename


class LocalProvider(BaseProvider):
    """
    Local genome provider.

    Give a local genome the genomepy treatment.
    """

    name = "Local"
    _cli_install_options = {
        "path_to_annotation": {
            "long": "path-to-annotation",
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

    def get_genome_download_link(self, path, mask=None, **kwargs):
        path = cleanpath(path)
        if not os.path.exists(path):
            raise FileNotFoundError(f"Local path to genome does not exist: {path}")
        return path

    def get_annotation_download_link(self, name: str, **kwargs) -> str:
        """
        Return a filepath to a matching annotation.

        Parameters
        ----------
        name : str
            genome name
        **kwargs: dict, optional:
            path_to_annotation : direct path to the gene annotation

        Returns
        -------
        str
            path

        Raises
        ------
        GenomeDownloadError
            if no functional path was found
        """
        path = kwargs.get("path_to_annotation")
        if path:
            path = cleanpath(path)
            if not os.path.exists(path):
                raise FileNotFoundError(
                    f"Local path to annotation does not exist: {path}"
                )
            ext = get_file_info(path)[0]
            if ext not in [".gtf", ".gff", ".gff3", ".bed"]:
                raise TypeError(
                    "Only (gzipped) gtf, gff and bed files are supported.\n"
                )
            return path

        paths = self.get_annotation_download_links(name)
        if paths:
            return paths[0]

        raise GenomeDownloadError(
            f"No gene annotations found for {get_genomename(name)}.\n"
        )

    def get_annotation_download_links(self, name, **kwargs):
        """Returns all files containing both name and an annotation extension"""
        name = cleanpath(name)
        genome_dir = os.path.dirname(name)
        search_list = os.listdir(genome_dir)
        search_name = get_genomename(name)

        hits = []
        for ext in ["gtf", "gff", "gff3"]:
            # .*? = non greedy filler. (\.gz)? = optional .gz
            expr = f"{search_name}.*?\.{ext}(\.gz)?"  # noqa: W605
            for line in search_list:
                hit = re.search(expr, line, flags=re.IGNORECASE)
                if hit:
                    hit = os.path.join(genome_dir, hit[0])
                    hits.append(hit)

        return hits
