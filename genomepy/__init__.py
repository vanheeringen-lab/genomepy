"""Search, download and use genome FASTA files."""
import genomepy.functions
__version__ = '0.2.1'
__author__ = "Simon van Heeringen"

list_available_providers = genomepy.functions.list_available_providers
list_available_genomes = genomepy.functions.list_available_genomes
list_installed_genomes = genomepy.functions.list_installed_genomes
search = genomepy.functions.search
install_genome = genomepy.functions.install_genome
genome = genomepy.functions.genome
