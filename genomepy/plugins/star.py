import os
import re
import subprocess as sp
from shutil import rmtree
from genomepy.plugin import Plugin
from genomepy.utils import mkdir_p, cmd_ok, run_index_cmd


class StarPlugin(Plugin):
    def after_genome_download(self, genome, force=False):
        if not cmd_ok("STAR"):
            return

        # Create index dir
        index_dir = genome.props["star"]["index_dir"]
        index_name = genome.props["star"]["index_name"]
        if force:
            # Start from scratch
            rmtree(index_dir, ignore_errors=True)
        mkdir_p(index_dir)

        if not os.path.exists(index_name):
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
            cmd = "STAR --runMode genomeGenerate --genomeFastaFiles {} --genomeDir {} --outFileNamePrefix {}".format(
                fname, index_dir, index_dir
            )
            run_index_cmd("star", cmd)

            # Rezip genome if it was bgzipped
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
                os.path.dirname(genome.filename), "index", "star"
            ),
            "index_name": os.path.join(
                os.path.dirname(genome.filename), "index", "star", "SA"
            ),
        }
        return props
