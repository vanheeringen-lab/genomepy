import os
import re
import subprocess as sp
from shutil import rmtree
from genomepy.plugin import Plugin
from genomepy.utils import mkdir_p, cmd_ok, run_index_cmd


class Hisat2Plugin(Plugin):
    def after_genome_download(self, genome, force=False):
        if not cmd_ok("hisat2-build"):
            return

        # Create index dir
        index_dir = genome.props["hisat2"]["index_dir"]
        index_name = genome.props["hisat2"]["index_name"]
        if force:
            # Start from scratch
            rmtree(index_dir, ignore_errors=True)
        mkdir_p(index_dir)

        if not any(fname.endswith(".ht2") for fname in os.listdir(index_dir)):
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
            cmd = "hisat2-build {} {}".format(fname, index_name)
            run_index_cmd("hisat2", cmd)

            if bgzip:
                ret = sp.check_call(["bgzip", fname])
                if ret != 0:
                    raise Exception(
                        "Error bgzipping genome {}. ".format(fname)
                        + "Is tabix installed?"
                    )

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
