import os.path
import re
import sys
from urllib.request import urlopen
from genomepy.plugin import Plugin


class BlacklistPlugin(Plugin):
    base_url = "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/"
    http_dict = {
        "ce10": base_url + "ce10-C.elegans/ce10-blacklist.bed.gz",
        "dm3": base_url + "dm3-D.melanogaster/dm3-blacklist.bed.gz",
        "hg38": base_url + "hg38-human/hg38.blacklist.bed.gz",
        "hg19": base_url + "hg38-human/hg19.blacklist.bed.gz",
        "mm9": base_url + "mm9-mouse/mm9-blacklist.bed.gz",
        "mm10": base_url + "mm10-mouse/mm10.blacklist.bed.gz",
    }

    def after_genome_download(self, genome, force=False):
        props = self.get_properties(genome)
        fname = props["blacklist"]
        if force and os.path.exists(fname):
            # Start from scratch
            os.remove(fname)

        if not os.path.exists(fname):
            link = self.http_dict.get(genome.name)
            if link is None:
                sys.stderr.write("No blacklist found for {}\n".format(genome.name))
                return
            try:
                sys.stderr.write("Downloading blacklist {}\n".format(link))
                response = urlopen(link)
                with open(fname, "wb") as bed:
                    bed.write(response.read())
            except Exception as e:
                sys.stderr.write(e)
                sys.stderr.write(
                    "Could not download blacklist file from {}".format(link)
                )

    def get_properties(self, genome):
        props = {
            "blacklist": re.sub(".fa(.gz)?$", ".blacklist.bed.gz", genome.filename)
        }
        return props
