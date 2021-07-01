import os.path
import re
from glob import glob

from pyfaidx import Fasta

from genomepy.files import glob_ext_files, read_readme
from genomepy.plugins import get_active_plugins
from genomepy.utils import get_genomes_dir, safe


class Genome(Fasta):
    """
    Get pyfaidx Fasta object of genome

    Also generates an index file of the genome

    Parameters
    ----------
    name : str
        Genome name

    genomes_dir : str
        Genome installation directory

    Returns
    -------
    pyfaidx.Fasta object
    """

    # import methods
    from genomepy.genome.sequences import get_random_sequences, track2fasta

    # lazy attributes (loaded when called)
    # listed here for code autocompletion
    sizes: dict = None
    gaps: dict = None

    def __init__(self, name, genomes_dir=None, *args, **kwargs):
        self.genomes_dir = get_genomes_dir(genomes_dir, check_exist=False)
        self.name = os.path.basename(re.sub(".fa(.gz)?$", "", safe(name)))
        self.filename = self._parse_filename(name)
        super(Genome, self).__init__(self.filename, *args, **kwargs)

        # file paths
        self.genome_file = self.filename
        self.genome_dir = os.path.dirname(self.filename)
        self.index_file = self.genome_file + ".fai"
        self.sizes_file = self._check_support_file("sizes")
        self.gaps_file = self._check_support_file("gaps")
        self.annotation_gtf_file = self._check_annotation_file("gtf")
        self.annotation_bed_file = self._check_annotation_file("bed")
        self.readme_file = os.path.join(self.genome_dir, "README.txt")

        # genome attributes
        metadata, _ = read_readme(self.readme_file)
        self.tax_id = metadata["tax_id"]
        self.assembly_accession = metadata["assembly_accession"]

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
        path_name = os.path.abspath(os.path.expanduser(name))
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
