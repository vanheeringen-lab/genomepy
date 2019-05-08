import re
import sys
from genomepy.plugin import Plugin
try:
    from urllib.request import urlopen
except:
    from urllib import urlopen

class BlacklistPlugin(Plugin):

    http_dict = {
        "ce10" : "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/ce10-C.elegans/ce10-blacklist.bed.gz",
        "dm3"  : "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/dm3-D.melanogaster/dm3-blacklist.bed.gz",
        "hg38" : "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz",
        "mm9"  : "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm9-mouse/mm9-blacklist.bed.gz",
        "mm10" : "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz",
    }

    def after_genome_download(self, genome):
        props = self.get_properties(genome)
        fname = props["blacklist"]
       
        link = self.http_dict.get(genome.name)
        if link is None: 
            sys.stderr.write("No blacklist found for {}\n".format(genome.name))
        try:
            sys.stderr.write("Downloading blacklist {}\n".format(link))
            response = urlopen(link)
            with open(fname, "wb") as bed:
                bed.write(response.read())
        except Exception as e:
            print(e)
            print("Could not download blacklist file from {}".format(link))
    
    def get_properties(self, genome):
        props = {
               "blacklist": genome.filename.replace(".fa", ".blacklist.bed.gz"),
               }
        return props
