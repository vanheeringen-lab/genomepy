from genomepy.plugin import Plugin

class SizesPlugin(Plugin):

    def after_genome_download(self, genome):
        props = self.get_properties(genome)
        fname = props["sizes"]
        with open(fname, "w") as f:
            for seqname in genome.keys():
                f.write("{}\t{}\n".format(seqname, len(genome[seqname])))

    def get_properties(self, genome):
        props = {
               "sizes":genome.filename + ".sizes",
               }
        return props
