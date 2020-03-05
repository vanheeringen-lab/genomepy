import os
from shutil import rmtree
from genomepy.plugin import Plugin
from genomepy.utils import mkdir_p, cmd_ok, run_index_cmd, bgunzip_and_name, bgrezip


class StarPlugin(Plugin):
    def after_genome_download(self, genome, threads=1, force=False):
        if not cmd_ok("STAR"):
            return

        # Create index dir
        index_dir = genome.props["star"]["index_dir"]
        index_name = genome.props["star"]["index_name"]
        if force:
            # Start from scratch
            rmtree(index_dir, ignore_errors=True)
        mkdir_p(index_dir)

        if not os.path.exists(index_name):
            # unzip genome if zipped and return up-to-date genome name
            bgzip, fname = bgunzip_and_name(genome)

            # Create index
            cmd = "STAR --runMode genomeGenerate --runThreadN {} --genomeFastaFiles {} --genomeDir {} --outFileNamePrefix {}".format(
                threads, fname, index_dir, index_dir
            )
            run_index_cmd("star", cmd)

            # re-zip genome if it was unzipped prior
            bgrezip(bgzip, fname)

    def get_properties(self, genome):
        props = {
            "index_dir": os.path.join(
                os.path.dirname(genome.filename), "index", "star"
            ),
            "index_name": os.path.join(
                os.path.dirname(genome.filename), "index", "star", "SA"
            ),
        }
        return props
