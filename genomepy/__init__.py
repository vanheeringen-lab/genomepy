"""Search, download and use genome FASTA files."""
import genomepy.__about__
import genomepy.functions
__version__ = genomepy.__about__.__version__
__author__ = genomepy.__about__.__author__

list_available_providers = genomepy.functions.list_available_providers
list_available_genomes = genomepy.functions.list_available_genomes
list_installed_genomes = genomepy.functions.list_installed_genomes
search = genomepy.functions.search
install_genome = genomepy.functions.install_genome
genome = genomepy.functions.genome
