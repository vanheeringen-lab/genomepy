"""Utility functions."""
import errno
import os
import norns
import re
import sys
import urllib.request
import subprocess as sp
import tarfile

from glob import glob
from pyfaidx import Fasta
from tempfile import TemporaryDirectory

config = norns.config("genomepy", default="cfg/default.yaml")


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


def get_genome_dir(genome_dir=None, check_exist=True):
    """import genome_dir if none is given, and check validity"""
    if not genome_dir:
        genome_dir = config.get("genome_dir", None)
    if not genome_dir:
        raise norns.exceptions.ConfigError("Please provide or configure a genome_dir")

    genome_dir = os.path.expanduser(genome_dir)
    if not os.path.exists(genome_dir) and check_exist:
        raise FileNotFoundError(f"Genome_dir {genome_dir} does not exist!")

    return genome_dir


def safe(name):
    """Replace spaces with undescores."""
    return name.strip().replace(" ", "_")


def get_localname(name, localname=None):
    """
    Returns localname if localname is not None, else;
      if name is an url (URL provider): Returns parsed filename from name
      else: returns name
    """
    if localname is None:
        try:
            urllib.request.urlopen(name)
        except (IOError, ValueError):
            return safe(name)
        else:
            # try to get the name from the url
            name = name[name.rfind("/") + 1 :]
            return safe(name[: name.find(".")])
    else:
        return safe(localname)


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


def sanitize_annotation(genome, gtf_file=None, sizes_file=None, out_dir=None):
    """
    Matches the toplevel sequence names in annotation.gtf to those in genome.fa.

    The fasta and gtf formats dictate the 1st field in the genome header and gtf should match.
    In some cases the genome.fa has multiple names per header, the 1st field not matching those in the gtf.
    If this occurs a conversion table can be made to rename the sequence names in either file.

    This script changes the names in the gtf to match those in the genome.fa, if needed.

    genome: class object
        Genome class object (contains genome name and file paths)

    sizes_file: str , optional
        Path to the genome sizes file

    gtf_file: str , optional
        Path to the (gzipped) annotation.gtf file

    out_dir: str , optional
        Location of the genome (and the README.txt)
    """
    # set default file paths if none are given
    if out_dir is None:
        out_dir = os.path.dirname(genome.filename)
    if gtf_file is None:
        gtf_file = os.path.join(out_dir, genome.name + ".annotation.gtf.gz")
    if sizes_file is None:
        sizes_file = os.path.join(out_dir, genome.name + ".fa.sizes")

    # unzip gtf if zipped and return up-to-date name
    if gtf_file.endswith(".gz"):
        sp.check_call(f"gunzip -f {gtf_file}", shell=True)
        gtf_file = gtf_file[:-3]

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
                sp.check_call(f"gzip -f {gtf_file}", shell=True)

                readme = os.path.join(out_dir, "README.txt")
                with open(readme, "a") as f:
                    f.write("genome and annotation have matching sequence names.\n")
                return

    # generate a gtf with matching scaffold/chromosome IDs
    sys.stderr.write(
        "\nGenome and annotation do not have matching sequence names! Creating matching annotation files...\n"
    )

    # unzip genome if zipped and return up-to-date genome name
    bgzip, genome_file = bgunzip_and_name(genome)

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
            # re-zip genome if unzipped
            bgrezip(bgzip, genome_file)

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

    # re-zip genome if unzipped
    bgrezip(bgzip, genome_file)

    # remove old files
    with TemporaryDirectory(dir=out_dir) as tmpdir:
        old_gtf_file = os.path.join(tmpdir, os.path.basename(gtf_file))
        bed_file = gtf_file.replace("gtf", "bed")
        sp.check_call(f"mv {gtf_file} {old_gtf_file}", shell=True)
        sp.check_call(f"rm {bed_file}.gz", shell=True)

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
    sp.check_call(f"gzip -f {gtf_file}", shell=True)

    readme = os.path.join(out_dir, "README.txt")
    with open(readme, "a") as f:
        f.write("corrected annotation files generated succesfully.\n")
    return
