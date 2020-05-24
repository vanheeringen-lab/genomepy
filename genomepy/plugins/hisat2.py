import os
import subprocess as sp
from shutil import rmtree
from genomepy.plugin import Plugin
from genomepy.utils import (
    mkdir_p,
    cmd_ok,
    run_index_cmd,
    bgzip_and_name,
    gunzip_and_name,
    gzip_and_name,
)


class Hisat2Plugin(Plugin):
    def after_genome_download(self, genome, threads=1, force=False):
        index_name = genome.plugin["hisat2"]["index_name"]
        if not cmd_ok("hisat2-build") or (
            os.path.exists(f"{index_name}.1.ht2") and not force
        ):
            return

        index_dir = genome.plugin["hisat2"]["index_dir"]
        rmtree(index_dir, ignore_errors=True)
        mkdir_p(index_dir)

        # gunzip genome if bgzipped and return up-to-date genome name
        fname, bgzip = gunzip_and_name(genome.filename)

        # index command
        cmd = f"hisat2-build -p {threads} {fname} {index_name}"

        # if an annotation is present, generate a splice-aware index
        gtf_file = genome.annotation_gtf_file
        if gtf_file:
            # gunzip if gzipped
            gtf_file, gzip_file = gunzip_and_name(gtf_file)

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
                f"python3 {splice_script} {gtf_file} > {splice_file}", shell=True
            )

            exon_script = hisat_path + "_extract_exons.py"
            exon_file = os.path.join(genome.genome_dir, "exon_sites.txt")
            sp.check_call(f"python3 {exon_script} {gtf_file} > {exon_file}", shell=True)

            # re-gzip annotation if gunzipped
            gzip_and_name(gtf_file, gzip_file)

            # update index command with annotation
            cmd += f" --ss {splice_file} --exon {exon_file}"
        else:
            print("\nCreating Hisat2 index without annotation file.")

        # Create index
        run_index_cmd("hisat2", cmd)

        # re-bgzip genome if gunzipped
        bgzip_and_name(fname, bgzip)

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
