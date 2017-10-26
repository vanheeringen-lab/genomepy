import os
import sys
import subprocess as sp

from genomepy.plugin import Plugin
from genomepy.utils import mkdir_p

class Bowtie2Plugin(Plugin):
    def after_genome_download(self, genome):
        try: 
            sp.check_call("bowtie2-build", stderr=sp.PIPE)
        except sp.CalledProcessError:
            # return code of 1 with no argument
            pass
        except:
            sys.stderr.write("bowtie2-build not found, skipping\n")
            return
        
        # Create index dir
        index_dir = genome.props["bowtie2"]["index_dir"]
        index_name =  genome.props["bowtie2"]["index_name"] 
        mkdir_p(index_dir)

        # Create index
        ret = sp.check_call("bowtie2-build {} {}".format(genome.filename, index_name), shell=True)
        if ret != 0:
            sys.stderr.write("bowtie2-build index return non-zero")
                
    def get_properties(self, genome):
        props = {
            "index_dir": os.path.join(os.path.dirname(genome.filename), "index", "bowtie2"),
            "index_name": os.path.join(os.path.dirname(genome.filename), "index", "bowtie2", genome.name),
            }
        return props
