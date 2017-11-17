import os
import sys
import subprocess as sp

from genomepy.plugin import Plugin
from genomepy.utils import mkdir_p, cmd_ok, run_index_cmd

class GmapPlugin(Plugin):
    def after_genome_download(self, genome):
        if not cmd_ok("gmap_build"):
            return

        # Create index dir
        index_dir = genome.props["gmap"]["index_dir"]
        index_name =  genome.props["gmap"]["index_name"] 
        mkdir_p(index_dir)

        # Create index
        cmd = "gmap_build -D {} -d {} {}".format(
                index_dir, genome.name, genome.filename)
        run_index_cmd("gmap", cmd)
         
    def get_properties(self, genome):
        props = {
            "index_dir": os.path.join(
                os.path.dirname(genome.filename), "index" , "gmap"
                ),
            "index_name": os.path.join(
                os.path.dirname(genome.filename), "index", "gmap", genome.name
                ),
            }
        return props
