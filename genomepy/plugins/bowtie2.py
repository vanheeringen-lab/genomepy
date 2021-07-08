import os

from genomepy.plugins import Plugin
from genomepy.utils import cmd_ok, mkdir_p, rm_rf, run_index_cmd


class Bowtie2Plugin(Plugin):
    def after_genome_download(self, genome, threads=1, force=False):
        if not cmd_ok("bowtie2-build"):
            return

        # Create index dir
        index_dir = genome.plugin["bowtie2"]["index_dir"]
        index_name = genome.plugin["bowtie2"]["index_name"]
        if force:
            # Start from scratch
            rm_rf(index_dir)
        mkdir_p(index_dir)

        if not any(fname.endswith(".bt2") for fname in os.listdir(index_dir)):
            # Create index
            cmd = f"bowtie2-build --threads {threads} {genome.filename} {index_name}"
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
