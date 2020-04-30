"""Utility functions."""
import os
import norns
import re
import sys
import urllib.request
import subprocess as sp
import tarfile
import gzip
import shutil
import time

from glob import glob
from norns import exceptions
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

    if not os.path.exists(readme):
        return metadata, lines

    # if the readme exists, overwrite all metadata fields found
    with open(readme) as f:
        for line in f.readlines():
            if ": " in line:
                vals = line.strip().split(": ")
                metadata[vals[0].strip()] = (": ".join(vals[1:])).strip()
            else:
                line = line.strip("\n").strip(" ")
                # blank lines are allowed, but only one in a row
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
    os.makedirs(path, exist_ok=True)


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
        raise exceptions.ConfigError("Please provide or configure a genomes_dir")

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


def gunzip_and_name(fname):
    """
    Gunzips the file if gzipped (also works on bgzipped files)
    Returns up-to-date filename and if it was gunzipped
    """
    if fname.endswith(".gz"):
        with gzip.open(fname, "rb") as f_in:
            with open(fname[:-3], "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.unlink(fname)
        return fname[:-3], True
    return fname, False


def gzip_and_name(fname, gzip_file=True):
    """
    Gzip file if requested
    Returns up to date filename
    """
    if gzip_file:
        with open(fname, "rb") as f_in:
            with gzip.open(fname + ".gz", "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        fname += ".gz"
    return fname


def bgzip_and_name(fname, bgzip=True):
    """
    Bgzip file if requested
    Returns up to date filename
    """
    if bgzip:
        ret = sp.check_call(["bgzip", fname])
        fname += ".gz"
        if ret != 0:
            raise Exception(f"Error bgzipping genome {fname}. Is tabix installed?")
    return fname


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


def match_contigs(gtf_file, sizes_file):
    """Check if genome and annotation have (properly) matching contig names"""
    # grab the first non-comment line and extract the first element (contig name)
    with open(gtf_file, "r") as gtf:
        for line in gtf:
            if not line.startswith("#"):
                gtf_id = line.split("\t")[0]
                break
        else:
            sys.stderr.write(f"No genes found in {gtf_file}. Skipping sanitizing.")
            return None  # cannot continue

    # check if contig name matches the first element in any of the genome headers
    with open(sizes_file, "r") as sizes:
        for line in sizes:
            fa_id = line.split("\t")[0]
            if fa_id == gtf_id:
                return True  # contig names match

    return False  # contig names do not match


def contig_pos(gtf_file, genome_file):
    """
    Determine which element in the fasta header contains the
    location identifiers used in the annotation.gtf

    returns the position of the matching element, -1 for no match.
    """
    # all headers in the file should be formatted the same. so take the first
    header = []
    with open(genome_file, "r") as fa:
        for line in fa:
            if line.startswith(">"):
                header = line.strip(">\n").split(" ")
                break

    # go over the whole gtf until a match is found
    element_pos = -1
    with open(gtf_file, "r") as gtf:
        for line in gtf:
            if not line.startswith("#"):
                loc_id = line.strip().split("\t")[0]
                try:
                    element_pos = header.index(loc_id)
                    break
                except ValueError:
                    continue

    return element_pos


def contig_conversion(genome_file, element_pos):
    """
    build a conversion table
    returns a list of duplicate contig names
    """
    conversion_table = {}
    duplicate_contigs = []
    with open(genome_file, "r") as fa:
        for line in fa:
            if line.startswith(">"):
                line = line.strip(">\n").split(" ")
                if line[element_pos] in conversion_table:
                    duplicate_contigs.append(line[element_pos])
                conversion_table[line[element_pos]] = line[0]
    return conversion_table, list(set(duplicate_contigs))


def sanitize_gtf(old_gtf_file, new_gtf_file, conversion_table):
    """
    create a new gtf file with renamed contigs based on the conversion table
    returns a list of contigs not present in the genome
    """
    missing_contigs = []
    with open(old_gtf_file, "r") as oldgtf, open(new_gtf_file, "w") as newgtf:
        for line in oldgtf:
            line = line.split("\t")
            try:
                line[0] = conversion_table[line[0]]
            except KeyError:
                missing_contigs.append(line[0])
            line = "\t".join(line)
            newgtf.write(line)
    return list(set(missing_contigs))


def sanitize_annotation(genome):
    """
    Matches the contig names in annotation.gtf to those in genome.fa.

    The fasta and gtf formats dictate the 1st field in the genome header and gtf should match.
    In some cases the genome.fa has multiple names per header, the 1st field not matching those in the gtf.
    If this occurs a conversion table can be made to rename the sequence names in either file.

    This script changes the names in the gtf to match those in the genome.fa, if needed.
    A bed format annotation file is made from the gtf afterward.
    """

    def cleanup(reason, error=False):
        """rezip gtf if needed, update readme and give an error message"""
        gzip_and_name(gtf_file, gzip_file)

        _metadata, _lines = read_readme(readme)
        _metadata["sanitized annotation"] = reason
        write_readme(readme, _metadata, _lines)

        if error:
            sys.stderr.write(
                "\nWARNING: Cannot correct annotation files automatically!\n"
                + "Leaving original version in place.\n"
            )

    readme = genome.readme_file
    gtf_file = genome.annotation_gtf_file
    gtf_file, gzip_file = gunzip_and_name(gtf_file)

    # check if sanitizing is needed
    match = match_contigs(gtf_file, genome.sizes_file)
    # clean up if sanitizing is not possible/required
    if match:
        cleanup("not required")
        return
    elif match is None:
        cleanup("not possible", error=True)
        return

    sys.stderr.write(
        "\nGenome and annotation do not have matching sequence names! "
        "Creating matching annotation files...\n"
    )
    # bgzip, genome_file = bgunzip_and_name(genome)
    genome_file, bgzip = gunzip_and_name(genome.filename)

    # try to find the (gtf) contig position in genome header
    element_pos = contig_pos(gtf_file, genome_file)
    # clean up if sanitizing is not possible
    if element_pos == -1:
        bgzip_and_name(genome_file, bgzip)
        cleanup("not possible", error=True)
        return

    # build a conversion table and check for duplicate contigs
    conversion_table, duplicate_contigs = contig_conversion(genome_file, element_pos)
    # re-zip genome if it was unzipped earlier
    bgzip_and_name(genome_file, bgzip)
    # clean up if sanitizing is not possible
    if duplicate_contigs:
        sys.stderr.write(
            "\nGenome contains duplicate contig names!\n"
            "The following contigs were found duplicate found: "
            f"{', '.join(duplicate_contigs)}.\n"
        )
        cleanup("not possible", error=True)
        return

    # all checks passed! Generate a new set of annotation files
    with TemporaryDirectory(dir=genome.genome_dir) as tmpdir:
        # generate corrected gtf file
        new_gtf_file = os.path.join(tmpdir, os.path.basename(gtf_file))
        missing_contigs = sanitize_gtf(gtf_file, new_gtf_file, conversion_table)

        # generate corrected bed file from the gtf
        cmd = "gtfToGenePred {0} /dev/stdout | genePredToBed /dev/stdin {1}"
        new_bed_file = new_gtf_file.replace("gtf", "bed")
        sp.check_call(cmd.format(new_gtf_file, new_bed_file), shell=True)

        # overwrite old files and gzip new files in needed
        bed_file = gtf_file.replace("gtf", "bed")
        for src, dst, gzip_file in zip(
            [new_gtf_file, new_bed_file],
            [gtf_file, bed_file],
            [genome.annotation_gtf_file, genome.annotation_bed_file],
        ):
            os.replace(src, dst)
            gzip_file = True if gzip_file.endswith(".gz") else False
            gzip_and_name(dst, gzip_file)

    metadata, lines = read_readme(readme)
    metadata["sanitized annotation"] = "yes"
    if missing_contigs:
        lines += [
            "",
            "WARNING: annotation contains contigs not present in the genome!",
            "The following contigs were not found:",
            f"{', '.join(missing_contigs)}.",
        ]
        sys.stderr.write(
            "\nWARNING: annotation contains contigs not present in the genome!\n"
            "The following contigs were not found: "
            f"{', '.join(missing_contigs)}.\n"
            "These have been kept as-is. Other contigs were corrected!\n"
        )
    write_readme(readme, metadata, lines)
