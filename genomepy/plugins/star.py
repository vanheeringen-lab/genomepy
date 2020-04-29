import os
import subprocess as sp
from shutil import rmtree
from genomepy.plugin import Plugin
from genomepy.utils import (
    mkdir_p,
    cmd_ok,
    run_index_cmd,
    gunzip_and_name,
    bgzip_and_name,
)


class StarPlugin(Plugin):
    def after_genome_download(self, genome, threads=1, force=False):
        if not cmd_ok("STAR"):
            return

        # Create index dir
        index_dir = genome.plugin["star"]["index_dir"]
        index_name = genome.plugin["star"]["index_name"]
        if force:
            # Start from scratch
            rmtree(index_dir, ignore_errors=True)
        mkdir_p(index_dir)

        if not os.path.exists(index_name):
            # unzip genome if zipped and return up-to-date genome name
            fname, bgzip = gunzip_and_name(genome.filename)
            # unzip annotation if found
            annot = fname[:-2] + "annotation.gtf.gz"
            rezip = False
            if os.path.exists(annot):
                sp.check_call(f"gunzip -f {annot}", shell=True)
                rezip = True
            annot = annot[:-3]

            # Create index
            cmd = (
                f"STAR --runMode genomeGenerate --runThreadN {threads} "
                + f"--genomeFastaFiles {fname} --genomeDir {index_dir} "
                + f"--outFileNamePrefix {index_dir}"
            )
            if os.path.exists(annot):
                cmd += f" --sjdbGTFfile {annot}"
            else:
                print("\nGenerating STAR index without annotation file.\n\n")
            run_index_cmd("star", cmd)

            # re-zip genome if it was unzipped prior
            bgzip_and_name(fname, bgzip)
            # re-zip annotation if it was unzipped prior
            if rezip:
                sp.check_call(f"gzip {annot}", shell=True)

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
