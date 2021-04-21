"""Utility functions."""
import os
import itertools
import norns
import re
import sys
import urllib.error
import urllib.request
import requests
import subprocess as sp
import tarfile
import gzip
import shutil
import socket
import time
from typing import Optional, Tuple
from ftplib import FTP, all_errors
from glob import glob
from norns import exceptions
from pyfaidx import Fasta
from tqdm.auto import tqdm
from tempfile import mkdtemp

from genomepy.exceptions import GenomeDownloadError

config = norns.config("genomepy", default="cfg/default.yaml")


def download_file(url, filename):
    """
    Helper method handling downloading large files from `url` to `filename`.
    Returns a pointer to `filename`.
    """

    def decorated_pbar(total):
        """Displays a progress bar with download speeds in MB/s."""
        return tqdm(
            desc="Download",
            unit_scale=True,
            unit_divisor=1024,
            total=total,
            unit="B",
        )

    def write_n_update_pbar(data):
        pbar.update(len(data))
        f.write(data)

    if url.startswith("ftp"):
        ftp, target = connect_ftp_link(url)
        file_size = ftp.size(target)
        with open(filename, "wb") as f:
            pbar = decorated_pbar(file_size)

            ftp.retrbinary(f"RETR {target}", write_n_update_pbar)
            ftp.quit()  # logout

    else:
        r = requests.get(url, stream=True)
        file_size = int(r.headers.get("Content-Length", 0))
        with open(filename, "wb") as f:
            pbar = decorated_pbar(file_size)

            for chunk in r.iter_content(chunk_size=1024):
                if chunk:  # filter out keep-alive new chunks
                    write_n_update_pbar(chunk)

    pbar.close()  # stop updating pbar
    return filename


def connect_ftp_link(link, timeout=None):
    """
    anonymous login to ftp
    accepts link in the form of ftp://ftp.name.domain/... and ftp.name.domain/...
    """
    link = link.replace("ftp://", "")
    host = link.split("/")[0]
    target = link.split(host)[1]

    try:
        ftp = FTP(host, timeout=timeout)
    except socket.gaierror:
        raise GenomeDownloadError(f"FTP host not found: {host}")

    ftp.login()
    return ftp, target


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
    seqs = [k for k in fa.keys() if bool(re.search(regex, k)) is not invert_match]
    if len(seqs) == 0:
        raise Exception("No sequences left after filtering!")

    outfa = outfa if outfa else infa
    _fa_to_file(fa, seqs, outfa)
    rm_rf(f"{infa}.fai")  # old index
    return Fasta(outfa)


def mkdir_p(path):
    """ 'mkdir -p' in Python """
    path = os.path.expanduser(path)
    os.makedirs(path, exist_ok=True)


def rm_rf(path):
    """ 'rm -rf' in Python """
    path = os.path.expanduser(path)
    if os.path.isfile(path):
        try:
            os.unlink(path)
        except OSError:  # in case of NTFS related issues
            pass
    elif os.path.isdir(path):
        shutil.rmtree(path, ignore_errors=True)


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
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)

    # show a spinner while the command is running
    spinner = itertools.cycle(["-", "\\", "|", "/"])
    while p.poll() is None:
        sys.stdout.write(f"\rCreating {name} index... {next(spinner)}")
        time.sleep(0.15)
        sys.stdout.flush()
    sys.stdout.write("\n")

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

    ext: str , optional
        Filename extension (default: fa).

    Returns
    -------
        File names.
    """
    fnames = glob(os.path.join(dirname, f"*.{ext}*"))
    return [f for f in fnames if f.endswith((ext, f"{ext}.gz"))]


def get_genomes_dir(genomes_dir: str = None, check_exist: Optional[bool] = True) -> str:
    """import genomes_dir if none is given, and check validity"""
    if not genomes_dir:
        # backwards compatibility for "genome_dir" (this fixes issue #87)
        genomes_dir = config.get("genomes_dir", config.get("genome_dir", None))
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
        name = name.split("/")[-1]  # remove path
        name = name.replace(".gz", "")  # remove .gz
        name = os.path.splitext(name)[0]  # remove .fa/.fna/.fasta etc
        name = safe(name)  # remove spaces
        # remove unwanted substrings from the name (ex: _genomes or .est_)
        unwanted = [
            "genome",
            "genomic",
            "sequence",
            "dna",
            "cds",
            "pep",
            "transcript",
            "EST",
            "toplevel",
            "primary",
            "assembly",
        ]
        spacers = "(-?_?\.?)"  # noqa: W605
        for substring in unwanted:
            name = re.sub(
                f"{spacers}{substring}(s?){spacers}", "", name, flags=re.IGNORECASE
            )
        return name


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


def is_number(term):
    """check if term is a number. Returns bool"""
    if isinstance(term, int) or term.isdigit():
        return True


def try_except_pass(errors, func, *args):
    """try to return FUNC with ARGS, pass on ERRORS"""
    try:
        return func(*args)
    except errors:
        pass


def retry(func, tries, *args):
    """
    Retry functions with potential connection errors.

    *args are passed as variables to func.
    """
    _try = 1
    while _try <= tries:
        try:
            answer = func(*args)
            return answer
        except (urllib.error.URLError, socket.timeout):
            time.sleep(1)
            _try += 1


def check_url(url, max_tries=1, timeout=15):
    """Check if URL works. Returns bool"""

    def _check_url(_url, _timeout):
        if _url.startswith("ftp"):
            ftp, target = connect_ftp_link(_url, timeout=_timeout)
            try:
                listing = ftp.nlst(target)
            except all_errors:
                listing = []
            ftp.quit()  # logout
            if listing:
                return True
        else:
            ret = urllib.request.urlopen(_url, timeout=_timeout)
            if ret.getcode() == 200:
                return True
        return False

    return retry(_check_url, max_tries, url, timeout)


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
