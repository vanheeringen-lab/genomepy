import os.path
from shutil import move
from tempfile import mkdtemp

from genomepy.files import extracted_file
from genomepy.plugins import Plugin
from genomepy.utils import cmd_ok, rm_rf, run_index_cmd


class GmapPlugin(Plugin):
    def after_genome_download(self, genome, threads=1, force=False):
        if not cmd_ok("gmap_build"):
            return

        # Create index dir
        index_dir = genome.plugin["gmap"]["index_dir"]
        if force:
            # Start from scratch
            rm_rf(index_dir)

        if not os.path.exists(index_dir):
            # unzip genome if zipped and return up-to-date genome name
            with extracted_file(genome.filename) as fname:
                # gmap outputs a folder named genome.name
                # its content is moved to index dir, consistent with other plugins
                tmp_dir = mkdtemp(dir=".")
                # Create index
                cmd = f"gmap_build -D {tmp_dir} -d {genome.name} {fname}"
                run_index_cmd("gmap", cmd)

                # Move files to index_dir
                src = os.path.join(tmp_dir, genome.name)
                move(src, index_dir)
                rm_rf(tmp_dir)

    def get_properties(self, genome):
        props = {
            "index_dir": os.path.join(
                os.path.dirname(genome.filename), "index", "gmap"
            ),
            "index_name": os.path.join(
                os.path.dirname(genome.filename), "index", "gmap", genome.name
            ),
        }
        return props
