import os
import sys
import subprocess as sp

from genomepy.base import Plugin

class GmapPlugin(Plugin):
    def after_genome_download(self, genome):
        try: 
            sp.check_call("gmap_build", stdout=sp.PIPE)
        except sp.CalledProcessError:
            # bwa gives return code of 1 with no argument
            pass
        except:
            sys.stderr.write("gmap_build not found, skipping\n")
            return
        
        # Create index dir
        index_dir = genome.props["gmap"]["index_dir"]
        index_name =  genome.props["gmap"]["index_name"] 
        if not os.path.exists(index_dir):
            os.mkdir(index_dir)

        # Create index
        ret = sp.check_call("gmap_build -D {} -d {} {}".format(index_dir, genome.name, genome.filename), shell=True)
        if ret != 0:
            sys.stderr.write("gmap_build index return non-zero")
                
    def get_properties(self, genome):
        props = {
            "index_dir": os.path.join(os.path.dirname(genome.filename), "gmap_index"),
            "index_name": os.path.join(os.path.dirname(genome.filename), "gmap_index", genome.name),
            }
        return props
