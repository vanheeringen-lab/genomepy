"""Utility functions."""
import errno
import os
import re
import sys
import urllib.request
import subprocess as sp

from pyfaidx import Fasta


def generate_gap_bed(fname, outname):
    """ Generate a BED file with gap locations.

    Parameters
    ----------
    fname : str
        Filename of input FASTA file.

    outname : str
        Filename of output BED file.
    """
    f = Fasta(fname)
    with open(outname, "w") as bed:
        for chrom in f.keys():
            for m in re.finditer(r"N+", f[chrom][:].seq):
                bed.write("{}\t{}\t{}\n".format(chrom, m.start(0), m.end(0)))


def generate_fa_sizes(fname, outname):
    """ Generate a fa.sizes file.

    Parameters
    ----------
    fname : str
        Filename of input FASTA file.

    outname : str
        Filename of output BED file.
    """
    f = Fasta(fname)
    with open(outname, "w") as sizes:
        for seqname, seq in f.items():
            sizes.write("{}\t{}\n".format(seqname, len(seq)))


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
                "{} already exists, set force to True to overwrite".format(outfa)
            )

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


def mkdir_p(path):
    """ 'mkdir -p' in Python """
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def cmd_ok(cmd):
    """Returns True if cmd can be run."""
    try:
        sp.check_call(cmd, stderr=sp.PIPE, stdout=sp.PIPE)
    except sp.CalledProcessError:
        # bwa gives return code of 1 with no argument
        pass
    except Exception:
        sys.stderr.write("{} not found, skipping\n".format(cmd))
        return False
    return True


def run_index_cmd(name, cmd):
    """Run command, show errors if the returncode is non-zero."""
    sys.stderr.write("Creating {} index...\n".format(name))
    # Create index
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        sys.stderr.write("Index for {} failed\n".format(name))
        sys.stderr.write(stdout.decode("utf8"))
        sys.stderr.write(stderr.decode("utf8"))


def get_localname(name, localname):
    """
    Returns localname if localname is not None, else;
      if name is an url (URL provider): Returns parsed filename from name
      else: returns name
    """
    if localname is None:
        try:
            urllib.request.urlopen(name)
        except (IOError, ValueError):
            return name.replace(" ", "_")
        else:
            # try to get the name from the url
            name = name[name.rfind("/") + 1 :]
            return name[: name.find(".")]
    else:
        return localname.replace(" ", "_")


def bgunzip_and_name(genome):
    """
    If the genome is bgzipped it needs to be unzipped first
    Returns up to date fname and bgzip status
    """
    fname = genome.filename
    bgzip = False
    if fname.endswith(".gz"):
        ret = sp.check_call(["gunzip", fname])
        if ret != 0:
            raise Exception("Error gunzipping genome {}".format(fname))
        fname = re.sub(".gz$", "", fname)
        bgzip = True
    return bgzip, fname


def bgrezip(bgzip, fname):
    """Rezip genome if it was unzipped by bgunzip"""
    if bgzip:
        ret = sp.check_call(["bgzip", fname])
        if ret != 0:
            raise Exception(
                "Error bgzipping genome {}. ".format(fname) + "Is tabix installed?"
            )
    return
