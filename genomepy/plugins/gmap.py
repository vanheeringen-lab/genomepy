import os
import re
import subprocess as sp

from genomepy.plugin import Plugin
from genomepy.utils import mkdir_p, cmd_ok, run_index_cmd


class GmapPlugin(Plugin):
    def after_genome_download(self, genome):
        if not cmd_ok("gmap_build"):
            return

        # Create index dir
        index_dir = genome.props["gmap"]["index_dir"]
        mkdir_p(index_dir)

        # If the genome is bgzipped it needs to be unzipped first
        fname = genome.filename
        bgzip = False
        if fname.endswith(".gz"):
            ret = sp.check_call(["gunzip", fname])
            if ret != 0:
                raise Exception("Error gunzipping genome {}".format(fname))
            fname = re.sub(".gz$", "", fname)
            bgzip = True

        # Create index
        cmd = "gmap_build -D {} -d {} {}".format(
            index_dir, genome.name, genome.filename)
        run_index_cmd("gmap", cmd)

        if bgzip:
            ret = sp.check_call(["bgzip", fname])
            if ret != 0:
                raise Exception(
                    "Error bgzipping genome {}. ".format(fname) +
                    "Is tabix installed?")

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
