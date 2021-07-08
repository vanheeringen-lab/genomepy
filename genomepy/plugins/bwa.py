import os

from genomepy.plugins import Plugin
from genomepy.utils import cmd_ok, mkdir_p, rm_rf, run_index_cmd


class BwaPlugin(Plugin):
    def after_genome_download(self, genome, threads=1, force=False):
        if not cmd_ok("bwa"):
            return

        # Create index dir
        index_dir = genome.plugin["bwa"]["index_dir"]
        index_name = genome.plugin["bwa"]["index_name"]
        if force:
            # Start from scratch
            rm_rf(index_dir)
        mkdir_p(index_dir)

        if not any(fname.endswith(".bwt") for fname in os.listdir(index_dir)):
            # Create index
            if not os.path.exists(index_name):
                os.symlink(genome.filename, index_name)
            cmd = f"bwa index {index_name}"
            run_index_cmd("bwa", cmd)

    def get_properties(self, genome):
        props = {
            "index_dir": os.path.join(os.path.dirname(genome.filename), "index", "bwa"),
            "index_name": os.path.join(
                os.path.dirname(genome.filename),
                "index",
                "bwa",
                f"{genome.name}.fa",
            ),
        }
        return props
