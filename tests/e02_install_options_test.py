import genomepy
import gzip
import os
import pytest
import shutil

from tempfile import mkdtemp
from time import sleep
from platform import system

skip = False
if not skip:

    def validate_gzipped_gtf(fname):
        assert os.path.exists(fname)
        with gzip.open(fname, "r") as f:
            for line in f:
                line = line.decode()
                if line.startswith("#"):
                    continue
                vals = line.split("\t")
                assert 9 == len(vals)
                int(vals[3]), int(vals[4])
                break

    def validate_gzipped_bed(fname):
        assert os.path.exists(fname)
        with gzip.open(fname, "r") as f:
            for line in f:
                line = line.decode()
                if line.startswith("#"):
                    continue
                vals = line.split("\t")
                assert 12 == len(vals)
                int(vals[1]), int(vals[2])
                break

    @pytest.fixture(scope="module", params=["no-overwrite", "overwrite"])
    def force(request):
        return request.param

    @pytest.fixture(scope="module", params=["original_name", "use_localname"])
    def localname(request):
        return request.param

    @pytest.fixture(scope="module", params=["unzipped", "bgzipped"])
    def bgzip(request):
        return request.param

    def test_install_genome_options(
        force, localname, bgzip, genome="ASM2732v1", provider="NCBI"
    ):
        """Test force, localname and bgzip"""
        tmp = mkdtemp()
        force = False if force == "no-overwrite" else True
        localname = None if localname == "original_name" else "My_localname"
        bgzip = False if bgzip == "unzipped" else True

        genomepy.install_genome(
            genome,
            provider,
            genomes_dir=tmp,
            localname=localname,
            bgzip=bgzip,
            force=force,
        )

        # force test
        ext = ".fa.gz" if bgzip else ".fa"
        name = genomepy.utils.get_localname(genome, localname)
        path = os.path.join(tmp, name, name + ext)

        t0 = os.path.getmtime(path)
        # OSX rounds down getmtime to the second
        if system() != "Linux":
            sleep(1)
        genomepy.install_genome(
            genome,
            provider,
            genomes_dir=tmp,
            localname=localname,
            bgzip=bgzip,
            force=force,
        )

        t1 = os.path.getmtime(path)
        assert t0 != t1 if force else t0 == t1

        shutil.rmtree(tmp)

    def test_install_annotation_options(
        force, localname, genome="ASM14646v1", provider="NCBI"
    ):
        """Test force and localname with annotations"""
        tmp = mkdtemp()
        force = False if force == "no-overwrite" else True
        localname = None if localname == "original_name" else "My_localname"

        # create dummy fasta to skip download_genome step
        name = genomepy.utils.get_localname(genome, localname)
        path = os.path.join(tmp, name, name + ".fa")
        os.mkdir(os.path.dirname(path))
        with open(path, "w") as f:
            f.write(">Chr1\nAAAACCCCTTTTGGGG\n")
        genomepy.install_genome(
            genome,
            provider,
            genomes_dir=tmp,
            localname=localname,
            only_annotation=True,
            skip_sanitizing=True,
            force=False,
        )

        gtf = os.path.join(tmp, name, name + ".annotation.gtf.gz")
        validate_gzipped_gtf(gtf)

        bed = os.path.join(tmp, name, name + ".annotation.bed.gz")
        validate_gzipped_bed(bed)

        # force test
        t0 = os.path.getmtime(gtf)
        # OSX rounds down getmtime to the second
        if system() != "Linux":
            sleep(1)
        genomepy.install_genome(
            genome,
            provider,
            genomes_dir=tmp,
            localname=localname,
            only_annotation=True,
            skip_sanitizing=True,
            force=force,
        )

        t1 = os.path.getmtime(gtf)
        assert t0 != t1 if force else t0 == t1

        shutil.rmtree(tmp)
