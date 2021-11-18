"""Utility functions with files"""
import gzip
import os
import re
import shutil
import subprocess as sp
import tarfile
from contextlib import contextmanager
from glob import glob
from tempfile import TemporaryDirectory, mkdtemp
from typing import Optional, Tuple, Union
from zipfile import ZipFile

from tqdm.auto import tqdm

from genomepy.utils import cleanpath, rm_rf


def parse_file(fpath, skip: Union[tuple, str] = "#"):
    """basic file parsing"""
    fpath = cleanpath(fpath)
    with open(fpath) as lines:
        for line in lines:
            line = line.strip()
            if line.startswith(skip):
                continue
            yield line


def read_readme(readme: str) -> Tuple[dict, list]:
    """
    Parse a readme file.

    Parameters
    ----------
    readme: str
        filename

    Returns
    -------
    tuple
        metadata : dict
            genome metadata

        lines: list
            non-metadata text (such as regex info)
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
    """
    Create a new readme file with supplied information.

    Parameters
    ----------
    readme: str
        filename
    metadata: dict
        dict to write
    lines: list, optional
        additional lines of text to append
    """
    with open(readme, "w") as f:
        for k, v in metadata.items():
            print(f"{k}: {v}", file=f)
        if lines:
            for line in lines:
                print(line, file=f)


def update_readme(readme: str, updated_metadata: dict = None, extra_lines: list = None):
    """
    Update a readme file with supplied information.

    Parameters
    ----------
    readme: str
        filename
    updated_metadata: dict
        dict with updated data
    extra_lines: list, optional
        additional lines of text to append
    """
    metadata, lines = read_readme(readme)
    if updated_metadata:
        metadata = {**metadata, **updated_metadata}
    if extra_lines:
        lines = lines + extra_lines
    write_readme(readme, metadata, lines)


@contextmanager
def extracted_file(fname: str):
    """Context manager to work with (b)gzipped file."""
    new_fname = extract_gzip(fname)
    gzipped = True
    if new_fname is None:
        new_fname = fname
        gzipped = False

    try:
        yield new_fname
    finally:
        if gzipped:
            try:
                bgzip_and_name(new_fname)
            except Exception:
                gzip_and_name(new_fname)


def extract_archive(
    fname: str, outfile: Optional[str] = None, concat: bool = False
) -> Union[str, None]:
    """
    Extract files from an archive.

    Archive may be a gzipped or bgzipped file (.gz), a zipped file (.zip)
    or a tarball (.tar.gz). Optionally, if multiple files are present
    they may be concatenated into one file.

    Parameters
    ----------
    fname: str
        path to file

    outfile : str, optional
        Name of the output file.

    concat : bool, optional
        If zip or tar.gz contains multiple files, concatenate them.

    Returns
    -------
    tuple
        fname, str
            up to date filename
    """
    if fname.endswith((".tgz", ".tar.gz")):
        return extract_tarball(fname, outfile=outfile)
    elif fname.endswith(".gz"):
        return extract_gzip(fname, outfile=outfile)
    elif fname.endswith(
        ".zip",
    ):
        return extract_zip(fname, outfile=outfile, concat=concat)


def extract_tarball(fname, outfile=None, concat=True) -> Union[str, None]:
    """Convert tar of multiple FASTAs to one file."""
    fnames = []
    # Extract files to temporary directory
    tmp_dir = mkdtemp(dir=os.path.dirname(outfile))
    with tarfile.open(fname) as tar:
        tar.extractall(path=tmp_dir)
    for root, _, files in os.walk(tmp_dir):
        fnames += [os.path.join(root, fname) for fname in files]

    if len(fnames) > 1 and not concat:
        raise ValueError("tarball contains multiple files, but concat not specified!")

    # Concatenate (also works woth one file)
    with open(outfile, "w") as out:
        for infile in fnames:
            for line in open(infile):
                out.write(line)

    rm_rf(tmp_dir)

    return outfile


def extract_gzip(fname: str, outfile: Optional[str] = None) -> Union[str, None]:
    """
    Gunzips the file if gzipped.

    Also works on bgzipped files.

    Parameters
    ----------
    fname: str
        path to file

    outfile : str, optional
        Name of the output file.

    Returns
    -------
    tuple
        fname, str
            up to date filename, or None if failed
    """
    if not outfile:
        outfile = fname[:-3]

    if fname.endswith(".gz"):
        with gzip.open(fname, "rb") as f_in:
            with open(outfile, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.unlink(fname)
        return outfile


def extract_zip(
    fname: str, outfile: Optional[str] = None, concat: bool = False
) -> Union[str, None]:
    """
    Unzips the file if zipped.

    Parameters
    ----------
    fname: str
        path to file

    outfile : str, optional
        Name of the output file.

    concat : bool, optional
        If zip or tar.gz contains multiple files, concatenate them.
    Returns
    -------
    tuple
        fname, str
            up to date filename
    """
    if not outfile:
        outfile = fname[:-4]

    with ZipFile(fname, "r") as fzip, TemporaryDirectory() as tmpdir:
        fzip.extractall(path=tmpdir)
        fnames = glob(f"{tmpdir}/*")

        if len(fnames) > 1 and not concat:
            raise ValueError(
                "Downloaded zip file contains multiple files, but concat is not specified!"
            )

        # Concatenate (also works with one file)
        with open(outfile, "w") as out:
            for infile in fnames:
                for line in open(infile):
                    out.write(line)

    os.unlink(fname)

    return outfile


def gzip_and_name(fname, gzip_file=True) -> str:
    """
    Gzip file if requested.

    Parameters
    ----------
    fname: str
        path to file
    gzip_file: bool, optional
        whether to gzip the file

    Returns
    -------
    str
        up to date filename
    """
    if gzip_file:
        with open(fname, "rb") as f_in:
            with gzip.open(fname + ".gz", "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.unlink(fname)
        fname += ".gz"
    return fname


def bgzip_and_name(fname, bgzip_file=True) -> str:
    """
    Bgzip file if requested.

    Parameters
    ----------
    fname: str
        path to file
    bgzip_file: bool, optional
        whether to gzip the file

    Returns
    -------
    str
        up to date filename
    """
    if bgzip_file:
        ret = sp.check_call(["bgzip", fname])
        fname += ".gz"
        if ret != 0:
            raise Exception(f"Error bgzipping genome {fname}. Is tabix installed?")
    return fname


def _open(fname: str, mode: Optional[str] = "r"):
    """
    Return a function to open a (gzipped) file.

    Parameters
    ----------
    fname: str
        (gzipped) file path
    mode: str
        (r)ead or (w)rite.

    Returns
    -------
    generator
        opened file
    """
    if mode not in ["r", "w"]:
        raise ValueError("mode must be either 'r' or 'w'.")

    if fname.endswith(".gz"):
        return gzip.open(fname, mode + "t")
    return open(fname, mode)


def get_file_info(fname) -> Tuple[str, bool]:
    """
    Returns the lower case file type of a file, and if it is (g)zipped

    Parameters
    ----------
    fname: str
        filename

    Returns
    -------
    tuple
        ext: str
            fname file type
        is_compressed: bool
            whether fname was (g)zipped
    """
    fname = fname.lower()
    is_compressed = False
    if fname.endswith((".tgz", ".tar.gz")):
        is_compressed = True
        fname = re.sub(r"\.(tgz|tar\.gz)$", "", fname)
    elif fname.endswith(".gz"):
        is_compressed = True
        fname = fname[:-3]
    elif fname.endswith(".zip"):
        is_compressed = True
        fname = fname[:-4]
    split = os.path.splitext(fname)
    return split[1], is_compressed


def glob_ext_files(dirname, ext="fa") -> list:
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
    list
        file names
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
    shutil.move(infa, old_fname)

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
    shutil.move(new_fname, outfa if outfa else infa)
    rm_rf(tmp_dir)

    return excluded_contigs


def filter_fasta(
    infa: str,
    outfa: str = None,
    regex: str = ".*",
    invert_match: Optional[bool] = False,
) -> list:
    """
    Filter fasta file based on regex.

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
    list
        removed contigs
    """
    pattern = re.compile(regex)

    def keep(header):
        return bool(pattern.search(header)) is not invert_match

    return _apply_fasta_regex_func(infa, keep, outfa)
