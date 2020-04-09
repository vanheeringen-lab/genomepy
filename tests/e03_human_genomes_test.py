import genomepy
import shutil
import pytest
import os

from tempfile import mkdtemp

travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"

skip = True
if not skip:

    @pytest.mark.skipif(travis, reason="slow")
    def test_ucsc_human():
        """Test UCSC.

        Download human genome from UCSC and retrieve a
        specific sequence.
        """
        tmp = mkdtemp()
        genomepy.install_genome("hg38", "UCSC", genome_dir=tmp)
        g = genomepy.Genome("hg38", genome_dir=tmp)
        seq = g["chr6"][166168664:166168679]
        assert str(seq) == "CCTCCTCGCTCTCTT"
        shutil.rmtree(tmp)

    @pytest.mark.skipif(travis, reason="slow")
    def test_ensembl_human():
        """Test Ensembl.

        Download human genome from Ensembl and retrieve a
        specific sequence.
        """
        tmp = mkdtemp()
        genomepy.install_genome("GRCh38.p13", "Ensembl", genome_dir=tmp)
        g = genomepy.Genome("GRCh38.p13", genome_dir=tmp)
        seq = g["6"][166168664:166168679]
        assert str(seq) == "CCTCCTCGCTCTCTT"
        shutil.rmtree(tmp)

    @pytest.mark.skipif(travis, reason="slow")
    def test_ncbi_human():
        """Test NCBI.

        Download human genome from NCBI and retrieve a
        specific sequence.
        """
        tmp = mkdtemp()
        genomepy.install_genome("GRCh38.p13", "NCBI", genome_dir=tmp)
        g = genomepy.Genome("GRCh38.p13", genome_dir=tmp)
        seq = g["6"][166168664:166168679]
        assert str(seq) == "CCTCCTCGCTCTCTT"
        shutil.rmtree(tmp)
