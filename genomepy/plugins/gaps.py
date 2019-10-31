import os.path
import re

from genomepy.plugin import Plugin


class GapsPlugin(Plugin):

    def after_genome_download(self, genome, force=False):
        props = self.get_properties(genome)
        fname = props["gaps"]

        if not os.path.exists(fname) or force is True:
            with open(fname, "w") as bed:
                for chrom in genome.keys():
                    for m in re.finditer(r'N+', genome[chrom][:].seq):
                        bed.write("{}\t{}\t{}\n".format(
                            chrom, m.start(0), m.end(0)))

    def get_properties(self, genome):
        props = {"gaps": genome.filename.replace(".fa", ".gaps.bed")}
        return props
