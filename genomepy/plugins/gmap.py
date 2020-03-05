import os.path
from shutil import move, rmtree
from tempfile import TemporaryDirectory
from genomepy.plugin import Plugin
from genomepy.utils import cmd_ok, run_index_cmd, bgunzip_and_name, bgrezip


class GmapPlugin(Plugin):
    def after_genome_download(self, genome, threads=1, force=False):
        if not cmd_ok("gmap_build"):
            return

        # Create index dir
        index_dir = genome.props["gmap"]["index_dir"]
        if force:
            # Start from scratch
            rmtree(index_dir, ignore_errors=True)

        if not os.path.exists(index_dir):
            # unzip genome if zipped and return up-to-date genome name
            bgzip, fname = bgunzip_and_name(genome)

            # gmap outputs a folder named genome.name
            # its content is moved to index dir, consistent with other plugins
            with TemporaryDirectory() as tmpdir:
                # Create index
                cmd = "gmap_build -D {} -d {} {}".format(tmpdir, genome.name, fname)
                run_index_cmd("gmap", cmd)

                # Move files to index_dir
                src = os.path.join(tmpdir, genome.name)
                move(src, index_dir)

            # re-zip genome if unzipped
            bgrezip(bgzip, fname)

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
