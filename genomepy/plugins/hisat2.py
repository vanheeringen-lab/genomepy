import os
import sys
import subprocess as sp

from genomepy.base import Plugin

class Hisat2Plugin(Plugin):
    def after_genome_download(self, genome):
        try: 
            sp.check_call("hisat2-build", stderr=sp.PIPE)
        except sp.CalledProcessError:
            # bwa gives return code of 1 with no argument
            pass
        except:
            sys.stderr.write("hisat2-build not found, skipping\n")
            return
        
        # Create index dir
        index_dir = genome.props["hisat2"]["index_dir"]
        index_name =  genome.props["hisat2"]["index_name"] 
        if not os.path.exists(index_dir):
            os.mkdir(index_dir)

        # Create index
        ret = sp.check_call("hisat2-build {} {}".format(genome.filename, index_name), shell=True)
        if ret != 0:
            sys.stderr.write("hisat2-build index return non-zero")
                
    def get_properties(self, genome):
        props = {
            "index_dir": os.path.join(os.path.dirname(genome.filename), "hisat2_index"),
            "index_name": os.path.join(os.path.dirname(genome.filename), "hisat2_index", genome.name),
            }
        return props
