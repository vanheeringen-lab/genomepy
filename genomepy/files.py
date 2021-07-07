import gzip
import os
import re
import shutil
import subprocess as sp
from glob import glob
from tempfile import mkdtemp
from typing import Optional, Tuple

from tqdm.auto import tqdm

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


def _apply_fasta_regex_func(infa, regex_func, outfa=None):
    """
    filter a Fasta using the regex function.

    infa: path to genome fasta

    regex_func: a function that takes a contig header and returns a bool

    outfa: path to output fasta. If None, infa is overwritten

    returns a list of excluded contigs
    """
    # move the original file to a tmp folder
    out_dir = os.path.dirname(infa)
    tmp_dir = mkdtemp(dir=out_dir)
    old_fname = os.path.join(tmp_dir, "original") if outfa is None else infa
    new_fname = os.path.join(tmp_dir, "filtered")
    os.rename(infa, old_fname)

    # perform the filtering
    excluded_contigs = []
    keep_contig = True
    with open(old_fname) as old, open(new_fname, "w") as new:
        for line in tqdm(old, desc="Filtering Fasta", unit_scale=1, unit=" lines"):
            if line[0] == ">":
                keep_contig = regex_func(line)
                if keep_contig is False:
                    excluded_contigs.append(line[1:].split(" ")[0].strip())
            if keep_contig:
                new.write(line)

    # move the filtered file to the original folder
    os.rename(new_fname, outfa if outfa else infa)
    rm_rf(tmp_dir)

    return excluded_contigs


def filter_fasta(
    infa: str,
    outfa: str = None,
    regex: str = ".*",
    invert_match: Optional[bool] = False,
) -> list:
    """Filter fasta file based on regex.

    Parameters
    ----------
    infa : str
        Filename of the input fasta file.

    outfa : str, optional
        Filename of the output fasta file. If None, infa is overwritten.

    regex : str, optional
        Regular expression used for selecting sequences.
        Matches everything if left blank.

    invert_match : bool, optional
        Select all sequence *not* matching regex if set.

    Returns
    -------
    a list of removed contigs
    """
    pattern = re.compile(regex)

    def keep(header):
        return bool(pattern.search(header)) is not invert_match

    return _apply_fasta_regex_func(infa, keep, outfa)
