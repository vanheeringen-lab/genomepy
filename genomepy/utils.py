"""Utility functions."""
import os
from pyfaidx import Fasta

def generate_sizes(name, genome_dir):
    """Generate a sizes file with length of sequences in FASTA file."""
    fa = os.path.join(genome_dir, name, "{}.fa".format(name))
    sizes = fa + ".sizes"
    g = Fasta(fa)
    with open(sizes, "w") as f:
        for seqname in g.keys():
            f.write("{}\t{}\n".format(seqname, len(g[seqname])))
