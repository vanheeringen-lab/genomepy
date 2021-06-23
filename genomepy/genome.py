import os.path
import re
from bisect import bisect
from glob import glob
from random import random

from pyfaidx import Fasta, Sequence

from genomepy.files import (
    generate_fa_sizes,
    generate_gap_bed,
    glob_ext_files,
    read_readme,
)
from genomepy.plugin import get_active_plugins
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

    def __init__(self, name, genomes_dir=None, *args, **kwargs):
        self.genomes_dir = get_genomes_dir(genomes_dir, check_exist=False)
        self.name = self._parse_name(name)
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
        self.sizes = {}
        self.gaps = {}
        metadata, _ = read_readme(self.readme_file)
        self.tax_id = metadata["tax_id"]
        self.assembly_accession = metadata["assembly_accession"]

    @property
    def sizes(self):
        return self.__sizes

    @sizes.setter
    def sizes(self, contig_sizes):
        self.__sizes = contig_sizes

    @sizes.getter
    def sizes(self):
        """Return sizes per contig when requested

        Returns
        -------
        contig_sizes : dict
            a dictionary with contigs as key and their lengths as values
        """
        if self.__sizes == {}:
            with open(self.sizes_file) as f:
                for line in f:
                    contig, length = line.strip().split("\t")
                    self.__sizes[contig] = int(length)
        return self.__sizes

    @property
    def gaps(self):
        return self.__gaps

    @gaps.setter
    def gaps(self, gap_sizes):
        self.__gaps = gap_sizes

    @gaps.getter
    def gaps(self):
        """Return gap sizes per chromosome when requested

        Returns
        -------
        gap_sizes : dict
            a dictionary with chromosomes as key and the total number of
            Ns as values
        """
        if self.__gaps == {}:
            with open(self.gaps_file) as f:
                for line in f:
                    chrom, start, end = line.strip().split("\t")
                    start, end = int(start), int(end)
                    self.__gaps[chrom] = self.__gaps.get(chrom, 0) + end - start
        return self.__gaps

    @property
    def plugin(self):
        """dict of all active plugins and their properties"""
        p = dict()
        for plugin in get_active_plugins():
            p[plugin.name()] = plugin.get_properties(self)
        return p

    @staticmethod
    def _parse_name(name: str) -> str:
        """extract a safe name from file path, url or regular names"""
        return os.path.basename(re.sub(".fa(.gz)?$", "", safe(name)))

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

    def _bed_to_seqs(self, track, stranded=False, extend_up=0, extend_down=0):
        bufsize = 10000
        with open(track) as fin:
            lines = fin.readlines(bufsize)
            for line in lines:
                if line.startswith("#") or line.startswith("track"):
                    continue

                vals = line.strip().split("\t")
                chrom, start, end = str(vals[0]), int(vals[1]), int(vals[2])
                name = f"{chrom}:{start}-{end}"

                # there might be more...
                starts = [start]
                ends = [end]

                # BED4: add name column to name
                if len(vals) >= 4:
                    name = " ".join((name, vals[3]))

                # BED5: check strandedness
                rc = False
                if stranded and len(vals) >= 6:
                    rc = vals[5] == "-"

                # BED12: get all blocks
                if len(vals) >= 12:
                    starts = [int(x) for x in vals[11].split(",")[:-1]]
                    sizes = [int(x) for x in vals[10].split(",")[:-1]]
                    starts = [start + x for x in starts]
                    ends = [start + size for start, size in zip(starts, sizes)]
                # convert to 1-based counting
                starts = [start + 1 for start in starts]

                # extend
                if extend_up:
                    if rc:
                        ends[-1] += extend_up
                    else:
                        starts[0] -= extend_up
                if extend_down:
                    if rc:
                        starts[0] -= extend_down
                    else:
                        ends[-1] += extend_down

                intervals = zip(starts, ends)
                seq = self.get_spliced_seq(chrom, intervals, rc)
                yield Sequence(name, seq.seq)

                # load more lines if needed
                lines += fin.readlines(1)

    def _region_to_seq(self, region, extend_up=0, extend_down=0):
        chrom, coords = region.strip().split(":")
        start, end = [int(c) for c in coords.split("-")]
        start += 1
        start -= extend_up
        end += extend_down
        seq = self.get_seq(chrom, start, end)
        return seq.seq

    def _regions_to_seqs(self, track, extend_up=0, extend_down=0):
        if isinstance(track, list):
            for region in track:
                name = region.strip()
                seq = self._region_to_seq(name, extend_up, extend_down)
                yield Sequence(name, seq)
        else:
            with open(track) as fin:
                bufsize = 10000
                lines = fin.readlines(bufsize)
                for region in lines:
                    name = region.strip()
                    seq = self._region_to_seq(name, extend_up, extend_down)
                    yield Sequence(name, seq)

                    # load more lines if needed
                    lines += fin.readlines()

    @staticmethod
    def get_track_type(track):
        # region_p example: "chr1:123-456"
        region_p = re.compile(r"^([^\s]+):(\d+)-(\d+)$")
        if not isinstance(track, (str, bytes)) and isinstance(track, (list, tuple)):
            if isinstance(track[0], (str, bytes)) and region_p.search(track[0]):
                return "interval"
        with open(track) as fin:
            line = fin.readline().strip()
        if region_p.search(line):
            return "interval"
        return "bed"

    def track2fasta(
        self, track, fastafile=None, stranded=False, extend_up=0, extend_down=0
    ):
        """
        Return a list of fasta sequences as Sequence objects
        as directed from the track(s).

        Params:
        track: list/region file/bed file
            region(s) you wish to translate to fasta.
            Example input files can be found in genomepy/tests/data/regions.*

        fastafile: bool , optional
            return Sequences as list or save to file? (default: list)

        stranded: bool , optional
            return sequences for both strands? Required BED6 (or higher) as input (default: False)

        extend up/down: int , optional
            extend the sequences to either side? (command is strand sensitive, default: 0)
        """
        track_type = self.get_track_type(track)
        if track_type == "interval":
            seqqer = self._regions_to_seqs(
                track, extend_up=extend_up, extend_down=extend_down
            )
        else:
            seqqer = self._bed_to_seqs(
                track, stranded=stranded, extend_up=extend_up, extend_down=extend_down
            )

        if fastafile:
            with open(fastafile, "w") as fout:
                for seq in seqqer:
                    fout.write(f"{seq.__repr__()}\n")
        else:
            return [seq for seq in seqqer]

    @staticmethod
    def _weighted_selection(list_of_tuples, n):
        """
        Selects n random elements from a list of (weight, item) tuples.
        Based on code snippet by Nick Johnson
        """
        cuml = []
        items = []
        total_weight = 0.0
        for weight, item in list_of_tuples:
            total_weight += weight
            cuml.append(total_weight)
            items.append(item)

        return [items[bisect(cuml, random() * total_weight)] for _ in range(n)]

    def get_random_sequences(
        self, n=10, length=200, chroms=None, max_n=0.1, outtype="list"
    ):
        """Return random genomic sequences.

        Parameters
        ----------
        n : int , optional
            Number of sequences to return.

        length : int , optional
            Length of sequences to return.

        chroms : list , optional
            Return sequences only from these chromosomes.

        max_n : float , optional
            Maximum fraction of Ns.

        outtype : string , optional
            return the output as list or string.
            Options: "list" or "string", default: "list".

        Returns
        -------
        coords : list of lists/strings
            List with [chrom, start, end] genomic coordinates.
            String with "chrom:start-end" genomic coordinates
            (can be used as input for track2fasta).
        """
        if not chroms:
            chroms = self.keys()

        # dict of chromosome sizes after subtracting the number of Ns
        sizes = dict(
            [(chrom, len(self[chrom]) - self.gaps.get(chrom, 0)) for chrom in chroms]
        )

        # list of (tuples with) chromosomes and their size
        # (if that size is long enough for random sequence selection)
        lengths = [
            (sizes[x], x)
            for x in chroms
            if sizes[x] / len(self[x]) > 0.1 and sizes[x] > 10 * length
        ]
        if len(lengths) == 0:
            raise Exception("No contigs of sufficient size were found.")

        # random list of chromosomes from lengths (can have duplicates)
        chroms = self._weighted_selection(lengths, n)

        coords = []
        retries = 100
        cutoff = length * max_n
        for chrom in chroms:
            for _ in range(retries):
                start = int(random() * (sizes[chrom] - length))
                end = start + length
                count_n = self[chrom][start:end].seq.upper().count("N")
                if count_n <= cutoff:
                    break
            else:
                raise Exception(
                    f"Random subset ran {retries} times, "
                    f"but could not find a sequence with less than {cutoff} N's in {chrom}.\n"
                    "You can specify contigs using the CHROMS argument."
                )

            # list output example ["chr1", 123, 456]
            coords.append([chrom, start, end])

        if outtype != "list":
            # bed output example: "chr1:123-456"
            for i, region in enumerate(coords):
                coords[i] = [f"{region[0]}:{region[1]}-{region[2]}"]

        return coords
