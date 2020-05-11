import os.path
import re
import sys
import zlib

from urllib.request import urlopen

from genomepy.plugin import Plugin


class BlacklistPlugin(Plugin):
    base_url = "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/"
    http_dict = {
        "ce10": base_url + "ce10-C.elegans/ce10-blacklist.bed.gz",
        "dm3": base_url + "dm3-D.melanogaster/dm3-blacklist.bed.gz",
        "hg38": base_url + "hg38-human/hg38.blacklist.bed.gz",
        "hg19": base_url
        + "hg19-human/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz",
        "mm9": base_url + "mm9-mouse/mm9-blacklist.bed.gz",
        "mm10": base_url + "mm10-mouse/mm10.blacklist.bed.gz",
        # for testing purposes
        "this was a triumph": "I'm making a note here: 'Huge success'",
    }

    def after_genome_download(self, genome, threads=1, force=False):
        fname = self.get_properties(genome)["blacklist"]
        if os.path.exists(fname) and not force:
            return

        link = self.http_dict.get(genome.name.split(".")[0])
        if link is None:
            sys.stderr.write(f"No blacklist found for {genome.name}\n")
            return

        sys.stderr.write(f"Downloading blacklist {link}\n")
        try:
            response = urlopen(link)
            with open(fname, "wb") as bed:
                # unzip the response with some zlib magic
                unzipped = zlib.decompress(response.read(), 16 + zlib.MAX_WBITS)
                bed.write(unzipped)
        except Exception as e:
            sys.stderr.write(str(e))
            sys.stderr.write(f"Could not download blacklist file from {link}")

    def get_properties(self, genome):
        props = {"blacklist": re.sub(".fa(.gz)?$", ".blacklist.bed", genome.filename)}
        return props
