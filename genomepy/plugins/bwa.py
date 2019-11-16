import os
from shutil import rmtree
from genomepy.plugin import Plugin
from genomepy.utils import mkdir_p, cmd_ok, run_index_cmd


class BwaPlugin(Plugin):
    def after_genome_download(self, genome, force=False):
        if not cmd_ok("bwa"):
            return

        # Create index dir
        index_dir = genome.props["bwa"]["index_dir"]
        index_name = genome.props["bwa"]["index_name"]
        if force:
            # Start from scratch
            rmtree(index_dir, ignore_errors=True)
        mkdir_p(index_dir)

        if not any(fname.endswith(".bwt") for fname in os.listdir(index_dir)):
            # Create index
            if not os.path.exists(index_name):
                os.symlink(genome.filename, index_name)

            cmd = "bwa index {}".format(index_name)
            run_index_cmd("bwa", cmd)

    def get_properties(self, genome):
        props = {
            "index_dir": os.path.join(os.path.dirname(genome.filename), "index", "bwa"),
            "index_name": os.path.join(
                os.path.dirname(genome.filename),
                "index",
                "bwa",
                "{}.fa".format(genome.name),
            ),
        }
        return props
