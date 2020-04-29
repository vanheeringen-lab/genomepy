import os.path
from shutil import move, rmtree
from tempfile import TemporaryDirectory
from genomepy.plugin import Plugin
from genomepy.utils import cmd_ok, run_index_cmd, gunzip_and_name, bgzip_and_name


class GmapPlugin(Plugin):
    def after_genome_download(self, genome, threads=1, force=False):
        if not cmd_ok("gmap_build"):
            return

        # Create index dir
        index_dir = genome.plugin["gmap"]["index_dir"]
        if force:
            # Start from scratch
            rmtree(index_dir, ignore_errors=True)

        if not os.path.exists(index_dir):
            # unzip genome if zipped and return up-to-date genome name
            fname, bgzip = gunzip_and_name(genome.filename)

            # gmap outputs a folder named genome.name
            # its content is moved to index dir, consistent with other plugins
            with TemporaryDirectory() as tmpdir:
                # Create index
                cmd = f"gmap_build -D {tmpdir} -d {genome.name} {fname}"
                run_index_cmd("gmap", cmd)

                # Move files to index_dir
                src = os.path.join(tmpdir, genome.name)
                move(src, index_dir)

            # re-zip genome if unzipped
            bgzip_and_name(fname, bgzip)

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
