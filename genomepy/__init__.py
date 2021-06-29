"""Search, download and use genome FASTA files."""
import sys

from loguru import logger

from genomepy.__about__ import __author__, __version__
from genomepy.annotation import Annotation
from genomepy.config import manage_config
from genomepy.exceptions import GenomeDownloadError
from genomepy.functions import (
    clean,
    install_genome,
    list_available_genomes,
    list_installed_genomes,
)
from genomepy.genome import Genome
from genomepy.plugins import manage_plugins
from genomepy.providers import list_providers
from genomepy.providers import search_all as search

# logger is a singleton, configuration here will be used module-wide
logger.remove()
logger.add(
    sys.stderr,
    format="<green>{time:HH:mm:ss}</green> <bold>|</bold> <blue>{level}</blue> <bold>|</bold> {message}",
    level="INFO",
)
