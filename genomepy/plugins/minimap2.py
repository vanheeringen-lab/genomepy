import os

from genomepy.plugins import Plugin
from genomepy.utils import cmd_ok, mkdir_p, rm_rf, run_index_cmd


class Minimap2Plugin(Plugin):
    def after_genome_download(self, genome, threads=1, force=False):
        if not cmd_ok("minimap2"):
            return

        # Create index dir
        index_dir = genome.plugin["minimap2"]["index_dir"]
        index_name = genome.plugin["minimap2"]["index_name"]
        if force:
            # Start from scratch
            rm_rf(index_dir)
        mkdir_p(index_dir)

        if not any(fname.endswith(".mmi") for fname in os.listdir(index_dir)):
            # Create index
            cmd = f"minimap2 -t {threads} -d {index_name} {genome.filename}"
            run_index_cmd("minimap2", cmd)

    def get_properties(self, genome):
        props = {
            "index_dir": os.path.join(
                os.path.dirname(genome.filename), "index", "minimap2"
            ),
            "index_name": os.path.join(
                os.path.dirname(genome.filename),
                "index",
                "minimap2",
                f"{genome.name}.mmi",
            ),
        }
        return props
