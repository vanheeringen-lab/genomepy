import os
from shutil import rmtree
from genomepy.plugin import Plugin
from genomepy.utils import mkdir_p, cmd_ok, run_index_cmd


class Bowtie2Plugin(Plugin):
    def after_genome_download(self, genome, force=False):
        if not cmd_ok("bowtie2-build"):
            return

        # Create index dir
        index_dir = genome.props["bowtie2"]["index_dir"]
        index_name = genome.props["bowtie2"]["index_name"]
        if force:
            # Start from scratch
            rmtree(index_dir, ignore_errors=True)
        mkdir_p(index_dir)

        if not any(fname.endswith(".bt2") for fname in os.listdir(index_dir)):
            # Create index
            cmd = "bowtie2-build {} {}".format(genome.filename, index_name)
            run_index_cmd("bowtie2", cmd)

    def get_properties(self, genome):
        props = {
            "index_dir": os.path.join(
                os.path.dirname(genome.filename), "index", "bowtie2"
            ),
            "index_name": os.path.join(
                os.path.dirname(genome.filename), "index", "bowtie2", genome.name
            ),
        }
        return props
