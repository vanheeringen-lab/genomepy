import os.path
import re
import sys

from bisect import bisect
from glob import glob
from pyfaidx import Fasta, Sequence
from random import random

from genomepy.plugin import get_active_plugins
from genomepy.provider import ProviderBase
from genomepy.utils import (
    read_readme,
    write_readme,
    get_genomes_dir,
    get_localname,
    glob_ext_files,
    generate_fa_sizes,
    generate_gap_bed,
    safe,
)


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

    def __init__(self, name, genomes_dir=None):
        self.genomes_dir = get_genomes_dir(genomes_dir, check_exist=False)
        self.name = self._parse_name(name)
        self.filename = self._parse_filename(name)
        super(Genome, self).__init__(self.filename)

        # file paths
        self.genome_file = self.filename
        self.genome_dir = os.path.dirname(self.filename)
        self.index_file = self.filename + ".fai"
        self.sizes_file = self.genome_file + ".sizes"
        self.gaps_file = os.path.join(self.genome_dir, self.name + ".gaps.bed")
        self.readme_file = os.path.join(self.genome_dir, "README.txt")

        # genome attributes
        self.sizes = {}  # populate using contig_sizes
        self.gaps = {}  # populate using gaps_sizes
        metadata = self._read_metadata()
        self.tax_id = metadata.get("tax_id")
        self.assembly_accession = metadata.get("assembly_accession")

    @property
    def sizes_file(self):
        return self.__sizes_file

    @sizes_file.setter
    def sizes_file(self, fname):
        if not os.path.exists(fname):
            generate_fa_sizes(self.genome_file, fname)
        self.__sizes_file = fname

    @property
    def gaps_file(self):
        return self.__gaps_file

    @gaps_file.setter
    def gaps_file(self, fname):
        if not os.path.exists(fname):
            generate_gap_bed(self.genome_file, fname)
        self.__gaps_file = fname

    @property
    def plugin(self):
        """dict of all active plugins and their properties"""
        p = dict()
        for plugin in get_active_plugins():
            p[plugin.name()] = plugin.get_properties(self)
        return p

    @property
    def annotation_gtf_file(self):
        return self.check_annotation_file("gtf")

    @property
    def annotation_bed_file(self):
        return self.check_annotation_file("bed")

    @staticmethod
    def _parse_name(name):
        """extract a safe name from file path, url or regular names"""
        return os.path.basename(re.sub(".fa(.gz)?$", "", get_localname(name)))

    def _parse_filename(self, name):
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

    def check_annotation_file(self, ext):
        """returns (gzipped) annotation file"""
        pattern = os.path.join(self.genome_dir, self.name + ".annotation.")
        file = glob(pattern + ext + "*")
        return file[0] if file else None

    @staticmethod
    def _update_provider(metadata):
        """check if provider is missing and try to update"""
        metadata["provider"] = "Unknown"
        url = metadata.get("genome url", "").lower()
        for provider in ["Ensembl", "UCSC", "NCBI"]:
            if provider.lower() in url:
                metadata["provider"] = provider
                break

    @staticmethod
    def _update_tax_id(metadata, provider=None, genome=None):
        """check if tax_id is missing and try to update"""
        taxid = "na"
        if genome:
            taxid = provider.genome_taxid(genome)
            taxid = str(taxid) if taxid != 0 else "na"
        metadata["tax_id"] = taxid

    @staticmethod
    def _update_assembly_accession(metadata, provider=None, genome=None):
        """check if assembly_accession is missing and try to update"""
        accession = "na"
        if genome:
            accession = provider.assembly_accession(genome)
        metadata["assembly_accession"] = accession

    def _update_metadata(self, metadata):
        """check if there is missing info that can be updated"""
        print(f"Updating metadata in README.txt", file=sys.stderr)
        if metadata.get("provider", "na") == "na":
            self._update_provider(metadata)

        known_provider = metadata["provider"] in ["Ensembl", "UCSC", "NCBI"]
        name = safe(metadata.get("original name", ""))
        missing_info = any(
            key not in metadata for key in ["tax_id", "assembly_accession"]
        )
        p = genome = None
        if known_provider and name and missing_info:
            p = ProviderBase.create(metadata["provider"])
            genome = p.genomes.get(name)

        if "tax_id" not in metadata:
            self._update_tax_id(metadata, p, genome)
        if "assembly_accession" not in metadata:
            self._update_assembly_accession(metadata, p, genome)

    def _read_metadata(self):
        """
        Read genome metadata from genome README.txt (if it exists).
        """
        metadata, lines = read_readme(self.readme_file)

        if (
            metadata.get("provider", "na") == "na"
            or "tax_id" not in metadata
            or "assembly_accession" not in metadata
        ):
            self._update_metadata(metadata)
            write_readme(self.readme_file, metadata, lines)

        return metadata

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
        region_p = re.compile(r"^(.+):(\d+)-(\d+)$")
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

    def contig_sizes(self):
        """Return sizes per contig.

        Returns
        -------
        contig_sizes : dict
            a dictionary with contigs as key and their lengths as values
        """
        with open(self.sizes_file) as f:
            for line in f:
                contig, length = line.strip().split("\t")
                self.sizes[contig] = length
        return self.sizes

    def gap_sizes(self):
        """Return gap sizes per chromosome.

        Returns
        -------
        gap_sizes : dict
            a dictionary with chromosomes as key and the total number of
            Ns as values
        """
        with open(self.gaps_file) as f:
            for line in f:
                chrom, start, end = line.strip().split("\t")
                start, end = int(start), int(end)
                self.gaps[chrom] = self.gaps.get(chrom, 0) + end - start
        return self.gaps

    @staticmethod
    def _weighted_selection(l, n):
        """
        Selects n random elements from a list of (weight, item) tuples.
        Based on code snippet by Nick Johnson
        """
        cuml = []
        items = []
        total_weight = 0.0
        for weight, item in l:
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
        if self.gaps is None:
            self.gap_sizes()

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
