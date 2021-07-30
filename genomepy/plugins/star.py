import os

from loguru import logger

from genomepy.files import extracted_file
from genomepy.plugins import Plugin
from genomepy.utils import cmd_ok, mkdir_p, rm_rf, run_index_cmd


class StarPlugin(Plugin):
    def after_genome_download(self, genome, threads=1, force=False):
        index_name = genome.plugin["star"]["index_name"]
        if not cmd_ok("STAR") or (os.path.exists(index_name) and not force):
            return

        index_dir = genome.plugin["star"]["index_dir"]
        rm_rf(index_dir)
        mkdir_p(index_dir)

        # gunzip genome if bgzipped and return up-to-date genome name
        with extracted_file(genome.filename) as fname:
            # index command
            cmd = (
                f"STAR --runMode genomeGenerate --runThreadN {threads} "
                + f"--genomeFastaFiles {fname} --genomeDir {index_dir} "
                + f"--outFileNamePrefix {index_dir}"
            )

            # if an annotation is present, generate a splice-aware index
            gtf_file = genome.annotation_gtf_file
            if gtf_file:
                with extracted_file(gtf_file) as _gtf_file:
                    # update index command with annotation
                    cmd += f" --sjdbGTFfile {_gtf_file}"

                    # Create index
                    run_index_cmd("star", cmd)
            else:
                logger.info("Creating STAR index without annotation file.")
                # Create index
                run_index_cmd("star", cmd)

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
