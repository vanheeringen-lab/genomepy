import gzip
import os
import re
import shutil
import subprocess as sp
import tarfile
from glob import glob
from tempfile import mkdtemp
from typing import Optional, Tuple

from pyfaidx import Fasta

from genomepy.utils import rm_rf


def read_readme(readme: str) -> Tuple[dict, list]:
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
        "sanitized annotation": "na",
        "genomepy version": "na",
        "date": "na",
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


def write_readme(readme: str, metadata: dict, lines: list = None):
    """Create a new readme with updated information"""
    with open(readme, "w") as f:
        for k, v in metadata.items():
            print(f"{k}: {v}", file=f)
        if lines:
            for line in lines:
                print(line, file=f)


def update_readme(readme: str, updated_metadata: dict = None, extra_lines: list = None):
    metadata, lines = read_readme(readme)
    if updated_metadata:
        metadata = {**metadata, **updated_metadata}
    if extra_lines:
        lines = lines + extra_lines
    write_readme(readme, metadata, lines)


def generate_gap_bed(fname, outname):
    """Generate a BED file with gap locations.

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
    """Generate a fa.sizes file.

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


def generate_annot(template, target, overwrite=False):
    """
    Create an annotation file type from the other file type.

    Parameters
    ----------
    template: str
        a GTF or BED filepath.

    target: str
        filepath to save the new annotation to.

    overwrite: bool, optional
        overwrite existing target file?
    """
    exts = os.path.basename(template.lower()).split(".")
    exts = [e for e in exts if e in ["gtf", "bed"]]
    if len(exts) == 0:
        raise ValueError("Template file must be in GTF or BED format.")
    template_ext = exts[-1]

    if not overwrite and os.path.exists(target):
        raise FileExistsError(f"{target} already exists! Set overwrite=True to ignore.")

    target_dir = os.path.dirname(target)
    tmp_dir = mkdtemp(dir=target_dir)
    tmp_target = os.path.join(tmp_dir, "new_annot")

    if template_ext == "bed":
        cmd = "bedToGenePred {0} /dev/stdout | genePredToGtf -source=genomepy file /dev/stdin {1}"
    else:
        cmd = "gtfToGenePred -ignoreGroupsWithoutExons {0} /dev/stdout | genePredToBed /dev/stdin {1}"

    # unzip template if needed
    template, is_unzipped = gunzip_and_name(template)
    # create new file
    sp.check_call(cmd.format(template, tmp_target), shell=True)
    # gzip if needed
    tmp_target = gzip_and_name(tmp_target, target.endswith(".gz"))
    _ = gzip_and_name(template, is_unzipped)

    os.replace(tmp_target, target)
    rm_rf(tmp_dir)


def tar_to_bigfile(fname, outfile):
    """Convert tar of multiple FASTAs to one file."""
    fnames = []
    # Extract files to temporary directory
    tmp_dir = mkdtemp(dir=os.path.dirname(outfile))
    with tarfile.open(fname) as tar:
        tar.extractall(path=tmp_dir)
    for root, _, files in os.walk(tmp_dir):
        fnames += [os.path.join(root, fname) for fname in files]

    # Concatenate
    with open(outfile, "w") as out:
        for infile in fnames:
            for line in open(infile):
                out.write(line)

    rm_rf(tmp_dir)


def gunzip_and_name(fname: str) -> (str, bool):
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
        os.unlink(fname)
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


def delete_extensions(directory: str, exts: list):
    """remove (gzipped) files in a directory matching any given extension"""
    for ext in exts:
        [rm_rf(f) for f in glob_ext_files(directory, ext)]


def _open(fname: str, mode: Optional[str] = "r"):
    """
    Return a function to open a (gzipped) file.

    fname: (gzipped) file path
    mode: (r)ead or (w)rite.
    """
    if mode not in ["r", "w"]:
        raise ValueError("mode must be either 'r' or 'w'.")

    if fname.endswith(".gz"):
        return gzip.open(fname, mode + "t")
    return open(fname, mode)


def file_len(fname):
    with _open(fname) as f:
        for i, _ in enumerate(f):  # noqa: B007
            pass
    return i + 1


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


def glob_ext_files(dirname, ext="fa"):
    """
    Return (gzipped) file names in directory containing the given extension.

    Parameters
    ----------
    dirname: str
        Directory name.

    ext: str , optional
        Filename extension (default: fa).

    Returns
    -------
        File names.
    """
    fnames = glob(os.path.join(dirname, f"*.{ext}*"))
    return [f for f in fnames if f.endswith((ext, f"{ext}.gz"))]


def is_genome_dir(dirname):
    """
    Check if a directory contains a fasta file of the same name

    Parameters
    ----------
    dirname : str
        Directory name

    Returns
    ------
    bool
    """
    genome_file = os.path.join(dirname, f"{os.path.basename(dirname)}.fa")
    return os.path.exists(genome_file) or os.path.exists(f"{genome_file}.gz")


def _fa_to_file(fasta: Fasta, contigs: list, filepath: str):
    tmp_dir = mkdtemp(dir=os.path.dirname(filepath))
    tmpfa = os.path.join(tmp_dir, "regex.fa")
    with open(tmpfa, "w") as out:
        for chrom in contigs:
            out.write(f">{fasta[chrom].name}\n")
            out.write(f"{fasta[chrom][:].seq}\n")
    os.rename(tmpfa, filepath)
    rm_rf(tmp_dir)


def filter_fasta(
    infa: str,
    regex: str = ".*",
    invert_match: Optional[bool] = False,
    outfa: str = None,
) -> Fasta:
    """Filter fasta file based on regex.

    Parameters
    ----------
    infa : str
        Filename of input fasta file.

    outfa : str, optional
        Filename of output fasta file.
        Overwrites infa if left blank.

    regex : str, optional
        Regular expression used for selecting sequences.
        Matches everything if left blank.

    invert_match : bool, optional
        Select all sequence *not* matching regex if set.

    Returns
    -------
        fasta : Fasta instance
            pyfaidx Fasta instance of newly created file
    """
    fa = Fasta(infa)
    pattern = re.compile(regex)
    seqs = [k for k in fa.keys() if bool(pattern.search(k)) is not invert_match]
    if len(seqs) == 0:
        raise Exception("No sequences left after filtering!")

    outfa = outfa if outfa else infa
    _fa_to_file(fa, seqs, outfa)
    rm_rf(f"{infa}.fai")  # old index
    return Fasta(outfa)
