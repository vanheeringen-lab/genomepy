import os.path
import re
import zlib
from urllib.request import urlopen

from loguru import logger

from genomepy.plugins import Plugin


class BlacklistPlugin(Plugin):
    # stanford_url = "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/"
    # encode_url = "https://www.encodeproject.org/files/"
    github_url = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/"
    http_dict = {
        "ce10": github_url + "ce10-blacklist.v2.bed.gz",
        "ce11": github_url + "ce11-blacklist.v2.bed.gz",
        "dm3": github_url + "dm3-blacklist.v2.bed.gz",
        "dm6": github_url + "dm6-blacklist.v2.bed.gz",
        "hg19": github_url + "hg19-blacklist.v2.bed.gz",
        "GRCh37": github_url + "hg19-blacklist.v2.bed.gz",
        "hg38": github_url + "hg38-blacklist.v2.bed.gz",
        "GRCh38": github_url + "hg38-blacklist.v2.bed.gz",
        "mm9": github_url + "Blacklist_v1/mm9-blacklist.bed.gz",  # no v2
        "mm10": github_url + "mm10-blacklist.v2.bed.gz",
        "GRCm38": github_url + "mm10-blacklist.v2.bed.gz",
        # for testing purposes
        "this was a triumph": "I'm making a note here: 'Huge success'",
    }

    def after_genome_download(self, genome, threads=1, force=False):
        fname = self.get_properties(genome)["blacklist"]
        if os.path.exists(fname) and not force:
            return

        link = self.http_dict.get(genome.name.split(".")[0])
        if link is None:
            logger.warning(f"No blacklist found for {genome.name}")
            return

        logger.info(f"Downloading blacklist {link}")
        try:
            response = urlopen(link)
            with open(fname, "wb") as bed:
                # unzip the response with some zlib magic
                unzipped = zlib.decompress(response.read(), 16 + zlib.MAX_WBITS)
                bed.write(unzipped)
        except Exception as e:
            logger.error(str(e))
            logger.error(f"Could not download blacklist file from {link}")

        # convert UCSC format to Ensembl/NCBI format
        if "1" in genome.sizes:
            os.rename(fname, fname + ".tmp")
            with open(fname + ".tmp") as old_bed, open(fname, "w") as new_bed:
                for line in old_bed:
                    line = line[3:]  # remove "chr" from the start of each line
                    new_bed.write(line)
            os.unlink(fname + ".tmp")

    def get_properties(self, genome):
        props = {"blacklist": re.sub(".fa(.gz)?$", ".blacklist.bed", genome.filename)}
        return props
