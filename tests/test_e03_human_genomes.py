from tempfile import mkdtemp

import pytest

import genomepy
import genomepy.utils
from tests import travis

skip = True
if not skip:

    @pytest.mark.skipif(travis, reason="slow")
    def test_ucsc_human():
        """Test UCSC.

        Download human genome from UCSC and retrieve a
        specific sequence.
        """
        tmp = mkdtemp()
        genomepy.install_genome("hg38", "UCSC", genomes_dir=tmp)
        g = genomepy.Genome("hg38", genomes_dir=tmp)
        seq = g["chr6"][166168664:166168679]
        assert str(seq) == "CCTCCTCGCTCTCTT"
        genomepy.utils.rm_rf(tmp)

    @pytest.mark.skipif(travis, reason="slow")
    def test_ensembl_human():
        """Test Ensembl.

        Download human genome from Ensembl and retrieve a
        specific sequence.
        """
        tmp = mkdtemp()
        genomepy.install_genome("GRCh38.p13", "Ensembl", genomes_dir=tmp)
        g = genomepy.Genome("GRCh38.p13", genomes_dir=tmp)
        seq = g["6"][166168664:166168679]
        assert str(seq) == "CCTCCTCGCTCTCTT"
        genomepy.utils.rm_rf(tmp)

    @pytest.mark.skipif(travis, reason="slow")
    def test_ncbi_human():
        """Test NCBI.

        Download human genome from NCBI and retrieve a
        specific sequence.
        """
        tmp = mkdtemp()
        genomepy.install_genome("GRCh38.p13", "NCBI", genomes_dir=tmp)
        g = genomepy.Genome("GRCh38.p13", genomes_dir=tmp)
        seq = g["6"][166168664:166168679]
        assert str(seq) == "CCTCCTCGCTCTCTT"
        genomepy.utils.rm_rf(tmp)
