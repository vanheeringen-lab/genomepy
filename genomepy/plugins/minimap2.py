import os
import sys
import subprocess as sp

from genomepy.base import Plugin

class Minimap2Plugin(Plugin):
    def after_genome_download(self, genome):
        try: 
            sp.check_call("minimap2", stderr=sp.PIPE)
        except sp.CalledProcessError:
            # bwa gives return code of 1 with no argument
            pass
        except:
            sys.stderr.write("minimap2 not found, skipping\n")
            return
        
        # Create index dir
        index_dir = genome.props["minimap2"]["index_dir"]
        index_name =  genome.props["minimap2"]["index_name"] 
        if not os.path.exists(index_dir):
            os.mkdir(index_dir)

        # Create index
        ret = sp.check_call("minimap2 -d {} {}".format(index_name, genome.filename), shell=True)
        if ret != 0:
            sys.stderr.write("minimap2 index return non-zero")
                
    def get_properties(self, genome):
        props = {
            "index_dir": os.path.join(os.path.dirname(genome.filename), "minimap2_index"),
            "index_name": os.path.join(os.path.dirname(genome.filename), "minimap2_index", "{}.mmi".format(genome.name)),
            }
        return props
