import os.path
from genomepy.plugin import Plugin


class SizesPlugin(Plugin):
    def after_genome_download(self, genome, force=False):
        props = self.get_properties(genome)
        fname = props["sizes"]

        if not os.path.exists(fname) or force:
            with open(fname, "w") as f:
                for seqname in genome.keys():
                    f.write("{}\t{}\n".format(seqname, len(genome[seqname])))

    def get_properties(self, genome):
        props = {"sizes": genome.filename.replace(".gz", "") + ".sizes"}
        return props
