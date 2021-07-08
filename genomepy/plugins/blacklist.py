import os.path
import re
import zlib
from urllib.request import urlopen

from loguru import logger

from genomepy.plugins import Plugin


class BlacklistPlugin(Plugin):
    stanford_url = "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/"
    encode_url = "https://www.encodeproject.org/files/"
    http_dict = {
        "ce10": stanford_url + "ce10-C.elegans/ce10-blacklist.bed.gz",
        "dm3": stanford_url + "dm3-D.melanogaster/dm3-blacklist.bed.gz",
        "hg19": stanford_url
        + "hg19-human/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz",
        "hg38": encode_url + "ENCFF356LFX/@@download/ENCFF356LFX.bed.gz",
        "GRCh38": encode_url + "ENCFF356LFX/@@download/ENCFF356LFX.bed.gz",
        "mm9": stanford_url + "mm9-mouse/mm9-blacklist.bed.gz",
        "mm10": stanford_url + "mm10-mouse/mm10.blacklist.bed.gz",
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
        if genome.name.split(".")[0] == "GRCh38":
            os.rename(fname, fname + ".tmp")
            with open(fname + ".tmp") as old_bed, open(fname, "w") as new_bed:
                for line in old_bed:
                    line = line[3:]  # remove "chr" from the start of each line
                    new_bed.write(line)
            os.unlink(fname + ".tmp")

    def get_properties(self, genome):
        props = {"blacklist": re.sub(".fa(.gz)?$", ".blacklist.bed", genome.filename)}
        return props
