import os
import subprocess as sp
from shutil import rmtree
from genomepy.plugin import Plugin
from genomepy.utils import mkdir_p, cmd_ok, run_index_cmd, bgunzip_and_name, bgrezip


class Hisat2Plugin(Plugin):
    def after_genome_download(self, genome, threads=1, force=False):
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
            # unzip genome if zipped and return up-to-date genome name
            bgzip, fname = bgunzip_and_name(genome)
            # unzip annotation if found
            annot = fname[:-2] + "annotation.gtf.gz"
            rezip = False
            if os.path.exists(annot):
                sp.check_call(f"gunzip -f {annot}", shell=True)
                rezip = True
            annot = annot[:-3]

            # generate splice and exon site files to enhance indexing
            splice_file = os.path.join(os.path.dirname(annot), "splice_sites.txt")
            exon_file = os.path.join(os.path.dirname(annot), "exon_sites.txt")
            if os.path.exists(annot):
                proc = sp.Popen("which hisat2", stdout=sp.PIPE, shell=True)
                splice_script = (
                    proc.stdout.read().decode("utf8").strip()
                    + "_extract_splice_sites.py"
                )
                sp.check_call(
                    f"python3 {splice_script} {annot} > {splice_file}", shell=True
                )

                exon_script = splice_script[:-24] + "_extract_exons.py"
                sp.check_call(
                    f"python3 {exon_script} {annot} > {exon_file}", shell=True
                )

            # Create index
            cmd = f"hisat2-build -p {threads} {fname} {index_name}"
            if os.path.exists(splice_file) and os.path.exists(exon_file):
                cmd += f" --ss {splice_file} --exon {exon_file}"
            else:
                print("\nGenerating Hisat2 index without annotation file.\n\n")
            run_index_cmd("hisat2", cmd)

            # re-zip genome if unzipped
            bgrezip(bgzip, fname)
            # re-zip annotation if it was unzipped prior
            if rezip:
                sp.check_call(f"gzip {annot}", shell=True)

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
