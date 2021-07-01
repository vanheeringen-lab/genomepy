import sys

from loguru import logger

from genomepy.__about__ import __author__, __version__
from genomepy.annotation import Annotation
from genomepy.caching import clean
from genomepy.config import manage_config
from genomepy.exceptions import GenomeDownloadError
from genomepy.functions import (
    install_genome,
    list_available_genomes,
    list_installed_genomes,
)
from genomepy.genome import Genome
from genomepy.plugins import manage_plugins
from genomepy.providers import Provider, list_online_providers, list_providers, search

# everything made available with `from genomepy import *`
__all__ = [
    "__author__",
    "__version__",
    "Provider",
    "Genome",
    "Annotation",
    "clean",
    "manage_config",
    "GenomeDownloadError",
    "install_genome",
    "list_available_genomes",
    "list_installed_genomes",
    "manage_plugins",
    "list_online_providers",
    "list_providers",
    "search",
]

# No traceback
# sys.tracebacklimit = 0

# logger is a singleton, configuration here will be used module-wide
logger.remove()
logger.add(
    sys.stderr,
    format="<green>{time:HH:mm:ss}</green> <bold>|</bold> <blue>{level}</blue> <bold>|</bold> {message}",
    level="INFO",
)
