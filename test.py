from genomepy.provider import ProviderBase
from genomepy import Genome
import sys


genome_name = "ce10"
genome = Genome(genome_name)
tax_id = genome.tax_id

#p = ProviderBase.create("UCSC")
#print(p.assembly_accession("ci3"))
#sys.exit()

p = ProviderBase.create("Ensembl")
name, accession, *rest = [row for row in p.search(tax_id)][0]
print(name, tax_id)
if accession == genome.assembly_accession:
    print(f"Ensembl {name} matches {genome_name} by accession")
else:
    print(f"Could not find a matching genome in Ensembl")
