import os
import sys
import subprocess as sp

from genomepy.base import Plugin

class BwaIndexPlugin(Plugin):
    def after_genome_download(self, genome):
        try: 
            sp.check_call("bwa", stderr=sp.PIPE)
        except sp.CalledProcessError:
            # bwa gives return code of 1 with no argument
            pass
        except:
            sys.stderr.write("bwa not found, skipping\n")
            return
        
        # Create index dir
        index_dir = genome.props["bwa_index"]["index_dir"]
        index_fa =  genome.props["bwa_index"]["index_name"] 
        if not os.path.exists(index_dir):
            os.mkdir(index_dir)
        if not os.path.exists(index_fa):
            os.symlink(genome.filename, index_fa)

        # Create index
        ret = sp.check_call("bwa index {}".format(index_fa), shell=True)
        if ret != 0:
            sys.stderr.write("bwa index return non-zero")
                
    def get_properties(self, genome):
        props = {
            "index_dir": os.path.join(os.path.dirname(genome.filename), "bwa_index"),
            "index_name": os.path.join(os.path.dirname(genome.filename), "bwa_index", "{}.fa".format(genome.name)),
            }
        return props
