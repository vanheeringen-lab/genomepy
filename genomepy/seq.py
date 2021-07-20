"""as_seqdict function & overloads"""
import os
import re
import sys
from functools import singledispatch
from io import TextIOWrapper
from tempfile import NamedTemporaryFile

try:
    import pybedtools
except ImportError:
    pass

import numpy as np
import pyfaidx
from Bio.SeqIO.FastaIO import SimpleFastaParser  # noqa: this is biopython

from genomepy.genome import Genome

# Regular expression to check for region (chr:start-end or genome@chr:start-end)
region_p = re.compile(r"^[^@]+@([^\s]+):(\d+)-(\d+)$")


@singledispatch
def as_seqdict(
    to_convert, genome=None, minsize=None  # noqa: arguments used dispaced functions
):
    """
    Convert input to a dictionary with name as key and sequence as value.

    If the input contains genomic coordinates, the genome needs to be
    specified. If minsize is specified all sequences will be checked if they
    are not shorter than minsize. If regions (or a region file) are used as
    the input, the genome can optionally be specified in the region using the
    following format: genome@chrom:start-end.

    Current supported input types include:
    * FASTA, BED and region files.
    * List or numpy.ndarray of regions.
    * pyfaidx.Fasta object.
    * pybedtools.BedTool object.

    Parameters
    ----------
    to_convert : list, str, pyfaidx.Fasta or pybedtools.BedTool
        Input to convert to FASTA-like dictionary
    genome : str, optional
        Genomepy genome name.
    minsize : int or None, optional
        If specified, check if all sequences have at least size minsize.

    Returns
    -------
    dict
        sequence names as key and sequences as value.
    """
    raise NotImplementedError(f"Not implement for {type(to_convert)}")


def _check_minsize(fa, minsize):
    """
    Raise ValueError if there is any sequence that is shorter than minsize.
    If minsize is None the size will not be checked.
    """
    if minsize is None:
        return fa

    for name, seq in fa.items():
        if len(seq) < minsize:
            raise ValueError(f"sequence {name} is shorter than {minsize}")

    return fa


def _genomepy_convert(to_convert, genome, minsize=None):
    """
    Convert a variety of inputs using track2fasta().
    """
    if genome is None:
        raise ValueError("input file is not a FASTA file, need a genome!")

    g = Genome(genome)
    tmpfile = NamedTemporaryFile()
    g.track2fasta(to_convert, tmpfile.name)

    fa = as_seqdict(tmpfile.name)
    return _check_minsize(fa, minsize)


def _as_seqdict_genome_regions(regions, minsize=None):
    """
    Accepts list of regions where the genome is encoded in the region,
    using the genome@chrom:start-end format.
    """
    genomic_regions = {}
    for region in regions:
        genome, region = region.split("@")
        genomic_regions.setdefault(genome, []).append(region)

    # test if all genomes are installed
    for genome in genomic_regions:
        Genome(genome)

    tmpfa = NamedTemporaryFile(mode="w", delete=False)
    for genome, g_regions in genomic_regions.items():
        g = Genome(genome)

        fa = g.track2fasta(g_regions)

        for seq in fa:
            seq.name = f"{genome}@{seq.name}"
            print(seq.__repr__(), file=tmpfa)

    tmpfa.flush()

    # Open tempfile and restore original sequence order
    fa = as_seqdict(tmpfa.name)
    fa = {region: fa[region] for region in regions}
    return _check_minsize(fa, minsize)


@as_seqdict.register(list)
def _as_seqdict_list(to_convert, genome=None, minsize=None):
    """
    Accepts list of regions as input.
    """
    if region_p.match(to_convert[0]):
        return _as_seqdict_genome_regions(to_convert, minsize)

    return _genomepy_convert(to_convert, genome, minsize)


@as_seqdict.register(TextIOWrapper)
def _as_seqdict_file_object(
    to_convert, genome=None, minsize=None  # noqa: arguments used dispaced functions
):
    """
    Accepts file object as input, should be a FASTA file.
    """
    fa = {x: y for x, y in SimpleFastaParser(to_convert)}
    return _check_minsize(fa, minsize)


@as_seqdict.register(str)
def _as_seqdict_filename(to_convert, genome=None, minsize=None):
    """
    Accepts filename as input.
    """
    if not os.path.exists(to_convert):
        raise ValueError("Assuming filename, but it does not exist")

    f = open(to_convert)
    fa = as_seqdict(f)

    if any(fa):
        return _check_minsize(fa, minsize)

    with open(to_convert) as f:
        while True:
            line = f.readline()
            if line == "":
                break
            if not line.startswith("#"):
                break

        if line == "":
            raise IOError(f"empty file {to_convert}")

        if region_p.match(line.strip()):
            regions = [region.strip() for region in [line] + f.readlines()]
            return _as_seqdict_genome_regions(regions, minsize=None)

    # Biopython parser resulted in empty dict
    # Assuming it's a BED or region file
    return _genomepy_convert(to_convert, genome, minsize)


@as_seqdict.register(pyfaidx.Fasta)
def _as_seqdict_pyfaidx(
    to_convert, genome=None, minsize=None  # noqa: arguments used dispaced functions
):
    """
    Accepts pyfaidx.Fasta object as input.
    """
    fa = {k: str(v) for k, v in to_convert.items()}
    return _check_minsize(fa, minsize)


@as_seqdict.register(np.ndarray)
def _as_seqdict_array(to_convert, genome=None, minsize=None):
    """
    Accepts numpy.ndarray with regions as input.
    """
    return as_seqdict(list(to_convert), genome, minsize)


# pybedtools is loaded
if "pybedtools" in sys.modules:
    if hasattr(pybedtools, "BedTool"):

        @as_seqdict.register(pybedtools.bedtool.BedTool)
        def _as_seqdict_bedtool(to_convert, genome=None, minsize=None):
            """
            Accepts pybedtools.BedTool as input.
            """
            return _genomepy_convert(
                ["{}:{}-{}".format(*f[:3]) for f in to_convert], genome, minsize
            )
