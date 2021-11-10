"""Genome class, modules & related functions"""
import os.path
import re
from glob import glob

from pyfaidx import Fasta

from genomepy.files import glob_ext_files, read_readme
from genomepy.genome.sequences import get_random_sequences as _get_random_sequences
from genomepy.genome.sequences import track2fasta as _track2fasta
from genomepy.plugins import get_active_plugins
from genomepy.utils import cleanpath, get_genomes_dir, safe

__all__ = ["Genome", "generate_fa_sizes", "generate_gap_bed"]


class Genome(Fasta):
    """
    pyfaidx Fasta object of a genome with additional attributes & methods.

    Generates a genome index file, sizes file and gaps file of the genome.

    Parameters
    ----------
    name : str
        Genome name

    genomes_dir : str, optional
        Genome installation directory

    Returns
    -------
    pyfaidx.Fasta
        An object that provides a pygr compatible interface.
    """

    # import methods
    get_random_sequences = _get_random_sequences
    track2fasta = _track2fasta

    # lazy attributes (loaded when called)
    # listed here for code autocompletion
    sizes: dict = None
    "contents of the sizes file: contigs and their lengths"
    gaps: dict = None
    "contents of the gaps file: contigs and the number of Ns contained"

    def __init__(self, name, genomes_dir=None, *args, **kwargs):
        self.name = safe(os.path.basename(re.sub(r"\.fa(\.gz)?$", "", name)))
        "genome name"
        self.genomes_dir = get_genomes_dir(genomes_dir, check_exist=False)
        "path to the genomepy genomes directory"
        self.filename = self._parse_filename(name)
        super(Genome, self).__init__(self.filename, *args, **kwargs)

        # file paths
        self.genome_file = self.filename
        "path to the genome fasta"
        self.genome_dir = os.path.dirname(self.filename)
        "path to the genome directory"
        self.index_file = self.genome_file + ".fai"
        "path to the genome index"
        self.sizes_file = self._check_support_file("sizes")
        "path to the chromosome sizes file"
        self.gaps_file = self._check_support_file("gaps")
        "path to the chromosome gaps file"
        self.annotation_gtf_file = self._check_annotation_file("gtf")
        "path to the gene annotation GTF file"
        self.annotation_bed_file = self._check_annotation_file("bed")
        "path to the gene annotation BED file"
        self.readme_file = os.path.join(self.genome_dir, "README.txt")
        "path to the README file"

        # genome attributes
        metadata, _ = read_readme(self.readme_file)
        self.tax_id = metadata["tax_id"]
        "genome taxonomy identifier"
        self.assembly_accession = metadata["assembly_accession"]
        "genome assembly accession"

    # lazy attributes
    def __getattribute__(self, name):
        val = super(Genome, self).__getattribute__(name)
        if val is not None:
            return val

        # if the attribute is None/empty, check if it is a lazy attribute
        if name == "sizes":
            val = {}
            with open(self.sizes_file) as f:
                for line in f:
                    contig, length = line.strip().split("\t")
                    val[contig] = int(length)
            setattr(self, name, val)

        elif name == "gaps":
            val = {}
            with open(self.gaps_file) as f:
                for line in f:
                    chrom, start, end = line.strip().split("\t")
                    val[chrom] = val.get(chrom, 0) + int(end) - int(start)
            setattr(self, name, val)

        return val

    # lazily update attributes if upstream attribute is updated
    def __setattr__(self, name, value):
        if name == "sizes_file":
            self.sizes = None  # noqa
        elif name == "gaps_file":
            self.gaps = None  # noqa
        super(Genome, self).__setattr__(name, value)

    @property
    def plugin(self):
        """dict of all active plugins and their properties"""
        p = dict()
        for plugin in get_active_plugins():
            p[plugin.name] = plugin.get_properties(self)
        return p

    def _parse_filename(self, name: str) -> str:
        """
        accepts path to a fasta file, path to a fasta folder, or
        the name of a genome (e.g. hg38).

        returns the abspath to the fasta file
        """
        path_name = cleanpath(name)
        if os.path.isfile(path_name):
            return path_name

        default_genome_dir = os.path.join(self.genomes_dir, self.name)
        for f in glob_ext_files(path_name) + glob_ext_files(default_genome_dir):
            if self.name + ".fa" in os.path.basename(f):
                return f

        raise FileNotFoundError(
            f"could not find {self.name}.fa(.gz) in genome_dir {default_genome_dir}"
        )

    def _check_support_file(self, ftype):
        """generate support file if missing/outdated"""
        ext = ".fa.sizes" if ftype == "sizes" else ".gaps.bed"
        func = generate_fa_sizes if ftype == "sizes" else generate_gap_bed

        fname = os.path.join(self.genome_dir, self.name + ext)
        if not os.path.exists(fname) or os.path.getmtime(fname) < os.path.getmtime(
            self.genome_file
        ):
            func(self.genome_file, fname)
        return fname

    def _check_annotation_file(self, ext):
        """returns (gzipped) annotation file"""
        pattern = os.path.join(self.genome_dir, self.name + ".annotation.")
        file = glob(pattern + ext + "*")
        return file[0] if file else None


def generate_gap_bed(fname, outname):
    """Generate a BED file with gap locations.

    Parameters
    ----------
    fname : str
        Filename of input FASTA file.

    outname : str
        Filename of output BED file.
    """
    f = Fasta(fname)
    with open(outname, "w") as bed:
        for chrom in f.keys():
            for m in re.finditer(r"N+", f[chrom][:].seq):
                bed.write(f"{chrom}\t{m.start(0)}\t{m.end(0)}\n")


def generate_fa_sizes(fname, outname):
    """Generate a fa.sizes file.

    Parameters
    ----------
    fname : str
        Filename of input FASTA file.

    outname : str
        Filename of output BED file.
    """
    f = Fasta(fname)
    with open(outname, "w") as sizes:
        for seqname, seq in f.items():
            sizes.write(f"{seqname}\t{len(seq)}\n")
