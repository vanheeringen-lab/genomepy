import os
from shutil import rmtree
from genomepy.plugin import Plugin
from genomepy.utils import mkdir_p, cmd_ok, run_index_cmd, bgunzip_and_name, bgrezip


class Hisat2Plugin(Plugin):
    def after_genome_download(self, genome, threads=1, force=False):
        if not cmd_ok("hisat2-build"):
            return

        # Create index dir
        index_dir = genome.props["hisat2"]["index_dir"]
        index_name = genome.props["hisat2"]["index_name"]
        if force:
            # Start from scratch
            rmtree(index_dir, ignore_errors=True)
        mkdir_p(index_dir)

        if not any(fname.endswith(".ht2") for fname in os.listdir(index_dir)):
            # unzip genome if zipped and return up-to-date genome name
            bgzip, fname = bgunzip_and_name(genome)

            # Create index
            cmd = "hisat2-build -p {} {} {}".format(threads, fname, index_name)
            run_index_cmd("hisat2", cmd)

            # re-zip genome if unzipped
            bgrezip(bgzip, fname)

    def get_properties(self, genome):
        props = {
            "index_dir": os.path.join(
                os.path.dirname(genome.filename), "index", "hisat2"
            ),
            "index_name": os.path.join(
                os.path.dirname(genome.filename), "index", "hisat2", genome.name
            ),
        }
        return props
