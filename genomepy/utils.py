"""Utility functions."""
import errno
import glob
import os
import re
import sys
import urllib.request
import subprocess as sp

from pyfaidx import Fasta
from tempfile import TemporaryDirectory


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


def glob_ext_files(dirname, ext="fa"):
    """
    Return (gzipped) file names in directory containing the given extension.

    Parameters
    ----------
    dirname: str
        Directory name.

    ext: str
        Filename extension (default: fa).

    Returns
    -------
        File names.
    """
    fnames = glob.glob(os.path.join(dirname, "*." + ext + "*"))
    return [fname for fname in fnames if fname.endswith(ext) or fname.endswith("gz")]


def sanitize_annotation(genome_dir, localname):
    """
    Matches the toplevel sequence names in annotation.gtf to those in the genome.fa.

    genome_dir:
        Save location for all genomes.

    localname:
        Custom name used for this genome.
    """
    out_dir = os.path.join(genome_dir, localname)

    genome_file = glob_ext_files(out_dir, "fa")[0]
    sizes_file = genome_file + ".sizes"
    gtf_file = os.path.join(out_dir, localname + ".annotation.gtf")

    sp.check_call("gunzip -f {}".format(gtf_file + ".gz"), shell=True)

    # Check if genome and annotation have matching chromosome/scaffold names
    with open(gtf_file, "r") as gtf:
        for line in gtf:
            if not line.startswith("#"):
                gtf_id = line.split("\t")[0]
                break

    with open(sizes_file, "r") as sizes:
        for line in sizes:
            fa_id = line.split("\t")[0]
            if fa_id == gtf_id:
                sp.check_call("gzip -f {}".format(gtf_file), shell=True)

                readme = os.path.join(out_dir, "README.txt")
                with open(readme, "a") as f:
                    f.write("genome and annotation have matching sequence names.\n")
                return

    # generate a gtf with matching scaffold/chromosome IDs
    sys.stderr.write(
        "Genome and annotation do not have matching sequence names! Creating matching annotation files...\n"
    )

    # determine which element in the fasta header contains the location identifiers used in the annotation.gtf
    header = []
    with open(genome_file, "r") as fa:
        for line in fa:
            if line.startswith(">"):
                header = line.strip(">\n").split(" ")
                break

    with open(gtf_file, "r") as gtf:
        for line in gtf:
            if not line.startswith("#"):
                loc_id = line.strip().split("\t")[0]
                try:
                    element = header.index(loc_id)
                    break
                except ValueError:
                    continue
        else:
            sys.stderr.write(
                "WARNING: Cannot correct annotation files automatically!\n"
                + "Leaving original version in place."
            )
            return

    # build a conversion table
    ids = {}
    with open(genome_file, "r") as fa:
        for line in fa:
            if line.startswith(">"):
                line = line.strip(">\n").split(" ")
                if line[element] not in ids.keys():
                    ids.update({line[element]: line[0]})

    # remove old files
    with TemporaryDirectory(dir=out_dir) as tmpdir:
        old_gtf_file = os.path.join(tmpdir, os.path.basename(gtf_file))
        bed_file = gtf_file.replace("gtf", "bed")
        sp.check_call("mv {} {}".format(gtf_file, old_gtf_file), shell=True)
        sp.check_call("rm {}".format(bed_file + ".gz"), shell=True)

        # rename the location identifier in the new gtf (using the conversion table)
        with open(old_gtf_file, "r") as oldgtf, open(gtf_file, "w") as newgtf:
            for line in oldgtf:
                line = line.split("\t")
                line[0] = ids[line[0]]
                line = "\t".join(line)
                newgtf.write(line)

    # generate a bed with matching scaffold/chromosome IDs
    cmd = (
        "gtfToGenePred {0} /dev/stdout | " "genePredToBed /dev/stdin {1} && gzip -f {1}"
    )
    sp.check_call(cmd.format(gtf_file, bed_file), shell=True)
    sp.check_call("gzip -f {}".format(gtf_file), shell=True)

    readme = os.path.join(out_dir, "README.txt")
    with open(readme, "a") as f:
        f.write("corrected annotation files generated succesfully.\n")
    return
