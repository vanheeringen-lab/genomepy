import os
import subprocess as sp

from loguru import logger

from genomepy.files import extracted_file
from genomepy.plugins import Plugin
from genomepy.utils import cmd_ok, mkdir_p, rm_rf, run_index_cmd


class Hisat2Plugin(Plugin):
    def after_genome_download(self, genome, threads=1, force=False):
        index_name = genome.plugin["hisat2"]["index_name"]
        if not cmd_ok("hisat2-build") or (
            os.path.exists(f"{index_name}.1.ht2") and not force
        ):
            return

        index_dir = genome.plugin["hisat2"]["index_dir"]
        rm_rf(index_dir)
        mkdir_p(index_dir)

        # gunzip genome if bgzipped and return up-to-date genome name
        with extracted_file(genome.filename) as fname:

            # index command
            cmd = f"hisat2-build -p {threads} {fname} {index_name}"

            # if an annotation is present, generate a splice-aware index
            gtf_file = genome.annotation_gtf_file
            print(gtf_file)
            if gtf_file:
                with extracted_file(gtf_file) as _gtf_file:
                    # generate splice and exon site files to enhance indexing
                    hisat_path = (
                        sp.Popen("which hisat2", stdout=sp.PIPE, shell=True)
                        .stdout.read()
                        .decode("utf8")
                        .strip()
                    )
                    splice_script = hisat_path + "_extract_splice_sites.py"
                    splice_file = os.path.join(genome.genome_dir, "splice_sites.txt")
                    sp.check_call(
                        f"python3 {splice_script} {_gtf_file} > {splice_file}",
                        shell=True,
                    )

                    exon_script = hisat_path + "_extract_exons.py"
                    exon_file = os.path.join(genome.genome_dir, "exon_sites.txt")
                    sp.check_call(
                        f"python3 {exon_script} {_gtf_file} > {exon_file}", shell=True
                    )

                    # update index command with annotation
                    cmd += f" --ss {splice_file} --exon {exon_file}"
                    # Create index
                    run_index_cmd("hisat2", cmd)
            else:
                logger.info("Creating Hisat2 index without annotation file.")

                # Create index
                run_index_cmd("hisat2", cmd)
            print(gtf_file)

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
