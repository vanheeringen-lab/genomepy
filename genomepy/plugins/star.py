import os
from shutil import rmtree
from genomepy.plugin import Plugin
from genomepy.utils import (
    mkdir_p,
    cmd_ok,
    run_index_cmd,
    gunzip_and_name,
    bgzip_and_name,
    gzip_and_name,
)


class StarPlugin(Plugin):
    def after_genome_download(self, genome, threads=1, force=False):
        index_name = genome.plugin["star"]["index_name"]
        if not cmd_ok("STAR") or (os.path.exists(index_name) and not force):
            return

        index_dir = genome.plugin["star"]["index_dir"]
        rmtree(index_dir, ignore_errors=True)
        mkdir_p(index_dir)

        # gunzip genome if bgzipped and return up-to-date genome name
        fname, bgzip = gunzip_and_name(genome.filename)

        # index command
        cmd = (
            f"STAR --runMode genomeGenerate --runThreadN {threads} "
            + f"--genomeFastaFiles {fname} --genomeDir {index_dir} "
            + f"--outFileNamePrefix {index_dir}"
        )

        # if an annotation is present, generate a splice-aware index
        gtf_file = genome.annotation_gtf_file
        gzip_file = False
        if gtf_file:
            # gunzip if gzipped
            gtf_file, gzip_file = gunzip_and_name(gtf_file)

            # update index command with annotation
            cmd += f" --sjdbGTFfile {gtf_file}"
        else:
            print("\nCreating STAR index without annotation file.")

        # Create index
        run_index_cmd("star", cmd)

        # re-bgzip genome if gunzipped
        bgzip_and_name(fname, bgzip)

        # re-gzip annotation if gunzipped
        if gtf_file:
            gzip_and_name(gtf_file, gzip_file)

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
