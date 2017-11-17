import os
import sys
import subprocess as sp

from genomepy.plugin import Plugin
from genomepy.utils import mkdir_p, cmd_ok, run_index_cmd

class BwaPlugin(Plugin):
    def after_genome_download(self, genome):
        if not cmd_ok("bwa"):
            return
        
        # Create index dir
        index_dir = genome.props["bwa"]["index_dir"]
        index_fa =  genome.props["bwa"]["index_name"] 
        mkdir_p(index_dir)

        if not os.path.exists(index_fa):
            os.symlink(genome.filename, index_fa)

        cmd = "bwa index {}".format(index_fa)
        run_index_cmd("bwa", cmd)
               
    def get_properties(self, genome):
        props = {
            "index_dir": os.path.join(
                os.path.dirname(genome.filename), "index", "bwa"
                ),
            "index_name": os.path.join(
                os.path.dirname(genome.filename), 
                "index","bwa", "{}.fa".format(genome.name)
                ),
            }
        return props
