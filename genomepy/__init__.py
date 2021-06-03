"""Search, download and use genome FASTA files."""
import sys

from loguru import logger

from genomepy.functions import (  # noqa: F401
    clean,
    install_genome,
    list_available_genomes,
    list_installed_genomes,
    manage_config,
    manage_plugins,
)
from genomepy.annotation import Annotation  # noqa: F401
from genomepy.genome import Genome  # noqa: F401
from genomepy.provider import Provider  # noqa: F401
from genomepy.exceptions import GenomeDownloadError  # noqa: F401
from genomepy.__about__ import __version__, __author__  # noqa: F401

list_available_providers = Provider.list_providers  # noqa: F401
search = Provider.search_all  # noqa: F401

# logger is a singleton, configuration here will be used module-wide
logger.remove()
logger.add(
    sys.stderr,
    format="<green>{time:HH:mm:ss}</green> <bold>|</bold> <blue>{level}</blue> <bold>|</bold> {message}",
    level="INFO",
)
