"""Utility functions."""
import errno
import os
import norns
import re
import sys
import urllib.request
import subprocess as sp
import tarfile
import time

from glob import glob
from pyfaidx import Fasta
from tempfile import TemporaryDirectory

config = norns.config("genomepy", default="cfg/default.yaml")


def read_readme(readme):
    """
    parse readme file

    Parameters
    ----------
    readme: str
        filename

    Returns
    -------
    tuple
        metadata : dict with genome metadata
        lines: list with non-metadata text (such as regex info)
    """
    metadata = {
        "name": "na",
        "provider": "na",
        "original name": "na",
        "original filename": "na",
        "assembly_accession": "na",
        "tax_id": "na",
        "mask": "na",
        "genome url": "na",
        "annotation url": "na",
        "sanitized annotation": "no",
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    lines = []

    # if the readme exists, overwrite all metadata fields found
    if os.path.exists(readme):
        with open(readme) as f:
            for line in f.readlines():
                if ": " in line:
                    vals = line.strip().split(": ")
                    metadata[vals[0].strip()] = (": ".join(vals[1:])).strip()
                else:
                    line = line.strip("\n").strip(" ")
                    if not (
                        line == ""
                        and len(lines) > 0
                        and lines[len(lines) - 1].strip() == ""
                    ):
                        lines.append(line)

    return metadata, lines


def write_readme(readme, metadata, lines):
    """Create a new readme with updated information"""
    with open(readme, "w") as f:
        for k, v in metadata.items():
            print(f"{k}: {v}", file=f)
        for line in lines:
            print(line, file=f)


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
                bed.write(f"{chrom}\t{m.start(0)}\t{m.end(0)}\n")


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
            sizes.write(f"{seqname}\t{len(seq)}\n")


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
            raise FileExistsError(
                f"{outfa} already exists, set force to True to overwrite"
            )

    filt_function = re.compile(regex).search
    fa = Fasta(infa, filt_function=filt_function)
    seqs = fa.keys()
    if v:
        original_fa = Fasta(infa)
        seqs = [s for s in original_fa.keys() if s not in seqs]
        fa = original_fa

    if len(seqs) == 0:
        raise Exception("No sequences left after filtering!")

    with open(outfa, "w") as out:
        for chrom in seqs:
            out.write(f">{fa[chrom].name}\n")
            out.write(f"{fa[chrom][:].seq}\n")

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
        sys.stderr.write(f"{cmd} not found, skipping\n")
        return False
    return True


def run_index_cmd(name, cmd):
    """Run command, show errors if the returncode is non-zero."""
    sys.stderr.write(f"Creating {name} index...\n")
    # Create index
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        sys.stderr.write(f"Index for {name} failed\n")
        sys.stderr.write(stdout.decode("utf8"))
        sys.stderr.write(stderr.decode("utf8"))


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
    fnames = glob(os.path.join(dirname, "*." + ext + "*"))
    return [
        fname for fname in fnames if fname.endswith(ext) or fname.endswith(ext + ".gz")
    ]


def get_genomes_dir(genomes_dir=None, check_exist=True):
    """import genomes_dir if none is given, and check validity"""
    if not genomes_dir:
        genomes_dir = config.get("genomes_dir", None)
    if not genomes_dir:
        raise norns.exceptions.ConfigError("Please provide or configure a genomes_dir")

    genomes_dir = os.path.abspath(os.path.expanduser(genomes_dir))
    if not os.path.exists(genomes_dir) and check_exist:
        raise FileNotFoundError(f"Genomes_dir {genomes_dir} does not exist!")

    return genomes_dir


def safe(name):
    """Replace spaces with undescores."""
    return name.strip().replace(" ", "_")


def get_localname(name, localname=None):
    """
    Returns the safe version of the given localname, if provided.
    If no localname is provided, return the safe version of the name.
    If the name is a working URL, return the safe version of the filename.
    """
    if localname:
        return safe(localname)
    try:
        urllib.request.urlopen(name)
    except (IOError, ValueError):
        return safe(name)
    else:
        # try to get the name from the url
        name = name[name.rfind("/") + 1 :]
        name = safe(name[: name.find(".fa")])
        # remove potential unwanted text from the name (ex: _genomes or .est_)
        unwanted = ["genome", "sequence", "cds", "pep", "transcript", "EST"]
        name = re.sub(
            r"(\.?_?){}(s?\.?_?)".format("(s?\.?_?)|".join(unwanted)),  # noqa: W605
            "",
            name,
            flags=re.IGNORECASE,
        )
        return name


def tar_to_bigfile(fname, outfile):
    """Convert tar of multiple FASTAs to one file."""
    fnames = []
    with TemporaryDirectory() as tmpdir:
        # Extract files to temporary directory
        with tarfile.open(fname) as tar:
            tar.extractall(path=tmpdir)
        for root, _, files in os.walk(tmpdir):
            fnames += [os.path.join(root, fname) for fname in files]

        # Concatenate
        with open(outfile, "w") as out:
            for infile in fnames:
                for line in open(infile):
                    out.write(line)
                os.unlink(infile)


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
            raise Exception(f"Error gunzipping genome {fname}")
        fname = re.sub(".gz$", "", fname)
        bgzip = True
    return bgzip, fname


def bgrezip(bgzip, fname):
    """Rezip genome if it was unzipped by bgunzip"""
    if bgzip:
        ret = sp.check_call(["bgzip", fname])
        if ret != 0:
            raise Exception(f"Error bgzipping genome {fname}. Is tabix installed?")
    return


def is_number(term):
    """check if term is a number. Returns bool"""
    if isinstance(term, int) or term.isdigit():
        return True


def check_url(url):
    """Check if URL works. Returns bool"""
    try:
        ret = urllib.request.urlopen(url)
        # check return code for http(s) urls
        if url.startswith("ftp") or ret.getcode() == 200:
            return True
    except urllib.request.URLError:
        return False


def read_url(url):
    """Read a text-based URL."""
    response = urllib.request.urlopen(url)
    data = response.read()
    text = data.decode("utf-8")
    return text


def get_file_info(fname):
    """
    Returns the lower case file type of a file, and if it is gzipped

    fname: str
        filename
    """
    fname = fname.lower()
    gz = False
    if fname.endswith(".gz"):
        gz = True
        fname = fname[:-3]
    split = os.path.splitext(fname)
    return split[1], gz
