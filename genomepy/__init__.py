import norns
import genomepy.functions
__version__ = '0.1.0'
__author__ = "Simon van Heeringen"

config = norns.config("genomepy", default="cfg/default.yaml")

list_available_genomes = genomepy.functions.list_available_genomes
list_installed_genomes = genomepy.functions.list_installed_genomes
search = genomepy.functions.search
install_genome = genomepy.functions.install_genome
genome = genomepy.functions.genome
