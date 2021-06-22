"""Search, download and use genome FASTA files."""
import sys

from loguru import logger

from genomepy.__about__ import __author__, __version__
from genomepy.annotation import Annotation
from genomepy.exceptions import GenomeDownloadError
from genomepy.functions import (
    clean,
    install_genome,
    list_available_genomes,
    list_installed_genomes,
    manage_config,
    manage_plugins,
)
from genomepy.genome import Genome
from genomepy.provider import Provider

list_available_providers = Provider.list_providers
search = Provider.search_all

# logger is a singleton, configuration here will be used module-wide
logger.remove()
logger.add(
    sys.stderr,
    format="<green>{time:HH:mm:ss}</green> <bold>|</bold> <blue>{level}</blue> <bold>|</bold> {message}",
    level="INFO",
)
