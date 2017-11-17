import os
import sys
import subprocess as sp

from genomepy.plugin import Plugin
from genomepy.utils import mkdir_p,  cmd_ok, run_index_cmd

class Hisat2Plugin(Plugin):
    def after_genome_download(self, genome):
        if not cmd_ok("hisat2-build"):
            return

        # Create index dir
        index_dir = genome.props["hisat2"]["index_dir"]
        index_name =  genome.props["hisat2"]["index_name"] 
        mkdir_p(index_dir)

        # Create index
        cmd = "hisat2-build {} {}".format(genome.filename, index_name)
        run_index_cmd("hisat2", cmd)
                
    def get_properties(self, genome):
        props = {
            "index_dir": os.path.join(
                os.path.dirname(genome.filename), "index", "hisat2"
                ),
            "index_name": os.path.join(
                os.path.dirname(genome.filename), "index", "hisat2", genome.name
                ),
            }
        return props
