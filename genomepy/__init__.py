"""
Genes and genomes at your fingertips!
"""

import sys

from loguru import logger

from genomepy.__about__ import __author__, __version__
from genomepy.annotation import Annotation, query_mygene
from genomepy.caching import clean
from genomepy.config import manage_config
from genomepy.functions import (
    head_annotations,
    install_genome,
    list_available_genomes,
    list_installed_genomes,
)
from genomepy.genome import Genome
from genomepy.plugins import manage_plugins
from genomepy.providers import Provider, list_online_providers, list_providers, search

# Public API objects
# everything made available with `from genomepy import *`
__all__ = [
    # modules
    "genome",
    "annotation",
    "providers",
    "files",
    "online",
    "utils",
    "exceptions",
    # objects
    "__author__",
    "__version__",
    "Genome",
    "Annotation",
    "Provider",
    "search",
    "head_annotations",
    "install_genome",
    "list_providers",
    "list_online_providers",
    "list_installed_genomes",
    "list_available_genomes",
    "manage_plugins",
    "manage_config",
    "clean",
    "query_mygene",
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
