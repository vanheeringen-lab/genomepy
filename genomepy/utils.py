"""Utility functions."""
import os
import re
from pyfaidx import Fasta

def generate_sizes(name, genome_dir):
    """Generate a sizes file with length of sequences in FASTA file."""
    fa = os.path.join(genome_dir, name, "{}.fa".format(name))
    sizes = fa + ".sizes"
    g = Fasta(fa)
    with open(sizes, "w") as f:
        for seqname in g.keys():
            f.write("{}\t{}\n".format(seqname, len(g[seqname])))

def filter_fasta(infa, outfa, regex=".*", v=False, force=False):
    """Filter fasta file based on regex.

    Parameters
    ----------
    infa : str
        Filename of input fasta file.
    
    outfa : str
        Filename of output fasta file. Cannot be the same as infa.

    regex : str, optional
        Regular expression used for selecting sequences.

    v : bool, optional
        If set to True, select all sequence *not* matching regex.

    force : bool, optional
        If set to True, overwrite outfa if it already exists.

    Returns
    -------
        fasta : Fasta instance
            pyfaidx Fasta instance of newly created file
    """
    if infa == outfa:
        raise ValueError("Input and output FASTA are the same file.")

    if os.path.exists(outfa):
        if force:
            os.unlink(outfa)
            if os.path.exists(outfa + ".fai"):
                os.unlink(outfa + ".fai")
        else:
            raise ValueError(
                    "{} already exists, set force to True to overwrite".format(outfa))
            
    filt_function = re.compile(regex).search
    fa = Fasta(infa, filt_function=filt_function)
    seqs = fa.keys()
    if v:
        original_fa = Fasta(infa)
        seqs = [s for s in original_fa.keys() if s not in seqs]
        fa = original_fa
    
    if len(seqs) == 0:
        raise ValueError("No sequences left after filtering!")

    with open(outfa, "w") as out:
        for chrom in seqs:
            out.write(">{}\n".format(fa[chrom].name))
            out.write("{}\n".format(fa[chrom][:].seq))

    return Fasta(outfa)

