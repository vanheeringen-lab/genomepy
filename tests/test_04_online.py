import os
import urllib.error
import urllib.request
from tempfile import NamedTemporaryFile

import pytest

import genomepy.online
from tests import linux, travis


def test_download_file():
    # HTTP
    tmp = NamedTemporaryFile().name
    assert not os.path.exists(tmp)
    url = "http://hgdownload.soe.ucsc.edu/goldenPath/ailMel1/bigZips/md5sum.txt"
    genomepy.online.download_file(url, tmp)
    assert os.path.exists(tmp)

    # FTP (doesn't work on Travis-Linux)
    if not (travis and linux):
        tmp = NamedTemporaryFile().name
        assert not os.path.exists(tmp)
        url = "ftp://ftp.ncbi.nlm.nih.gov//genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/README.txt"
        genomepy.online.download_file(url, tmp)
        assert os.path.exists(tmp)


def test_connect_ftp_link():
    if not (travis and linux):  # (doesn't work on Travis-Linux)
        # good FTP host
        ftp_link = "ftp://ftp.ncbi.nlm.nih.gov/genomes/README.txt"
        ftp, target = genomepy.online.connect_ftp_link(ftp_link)
        assert target == "genomes/README.txt"
        result = ftp.nlst(target)
        ftp.quit()  # logout
        assert result == [target]

        # bad FTP host
        with pytest.raises(genomepy.exceptions.GenomeDownloadError):
            genomepy.online.connect_ftp_link("ftp://not.an.ftp/at/all")


def test_read_url(
    url="http://ftp.xenbase.org/pub/Genomics/JGI/README", expected="The data"
):
    text = genomepy.online.read_url(url)
    assert text.startswith(expected)


def test_retry(capsys):
    # runs passed function
    txt = "hello world"
    genomepy.online.retry(print, 1, txt)
    captured = capsys.readouterr().out.strip()
    assert captured == txt

    # handles URLErrors
    def _offline_func():
        raise urllib.error.URLError("this function is offline")

    assert genomepy.online.retry(_offline_func, 1) is None


def test_check_url():
    # good URL
    assert genomepy.online.check_url(
        "https://mbdata.science.ru.nl/siebrenf/GRCh37/GRCh37/README.txt"
    )
    assert genomepy.online.check_url("http://ftp.xenbase.org/pub/Genomics/JGI/README")

    # bad URL:
    assert not genomepy.online.check_url("https://www.thiswebsiteisoffline.nl/")

    if not (travis and linux):  # (doesn't work on Travis-Linux)
        # good FTP
        ftp_link = (
            "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/"
            "GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.gff.gz"
        )
        assert genomepy.online.check_url(ftp_link)

        # bad FTP target
        assert not genomepy.online.check_url("ftp://ftp.ncbi.nlm.nih.gov/bad_target")
