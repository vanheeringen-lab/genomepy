"""Utility functions with internet connections"""
import socket
import urllib.error
from ftplib import FTP, all_errors, error_temp
from time import sleep
from typing import Tuple
from urllib.request import urlopen

import requests
from tqdm.auto import tqdm

from genomepy.exceptions import GenomeDownloadError


def download_file(url, filename) -> str:
    """
    Helper method handling downloading large files from `url` to `filename`.

    Parameters
    ----------
    url : str
        download target url
    filename : str
        file to download to

    Returns
    -------
    str
        filename
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


def connect_ftp_link(link, timeout=None) -> Tuple[FTP, str]:
    """
    Anonymous login to ftp.

    Accepts link in the form of ftp://ftp.name.domain/...
    and ftp.name.domain/...

    Parameters
    ----------
    link : str
        FTP link
    timeout : int, optional
        number of idle seconds before the connection closes

    Returns
    -------
    tuple
        ftp: FTP
            object with connection established
        target: str
            target file
    """
    link = link.replace("ftp://", "")
    host, target = link.split("/", 1)
    try:
        ftp = FTP(host, timeout=timeout)
    except socket.gaierror:
        raise GenomeDownloadError(f"FTP host not found: {host}")
    ftp.login()
    return ftp, target


def read_url(url):
    """Read a text-based URL."""
    response = urlopen(url)
    data = response.read()
    text = data.decode("utf-8")
    return text


def retry(func, tries, *args, **kwargs):
    """
    Retry functions with potential connection errors.

    args and kwargs are passed to func.
    """
    for n in range(tries):
        try:
            answer = func(*args, **kwargs)
            return answer
        except (urllib.error.HTTPError, error_temp):
            # HTTPError 404: URL not found
            # FTP error_temp 450: file not found
            return
        except all_errors + (urllib.error.URLError, socket.timeout):
            # connection errors: wait and try again
            if n < tries - 1:
                sleep(1)


def check_url(url, max_tries=1, timeout=15) -> bool:
    """Check if URL works. Returns bool"""

    def _check_url(_url, _timeout):
        if _url.startswith("ftp"):
            ftp, target = connect_ftp_link(_url, timeout=_timeout)
            listing = retry(ftp.nlst, 1, target)
            ftp.quit()  # logout
            if listing:
                return True
        else:
            ret = urlopen(_url, timeout=_timeout)
            if ret.getcode() == 200:
                return True

        return False

    return retry(_check_url, max_tries, url, timeout)
