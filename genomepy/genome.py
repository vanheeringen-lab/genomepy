import os.path
import re
import sys

from bisect import bisect
from pyfaidx import Fasta, Sequence
from random import random

from genomepy.plugin import get_active_plugins
from genomepy.provider import ProviderBase
from genomepy.utils import (
    read_readme,
    write_readme,
    get_genome_dir,
    get_localname,
    glob_ext_files,
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

    genome_dir : str
        Genome installation directory

    Returns
    -------
    pyfaidx.Fasta object
    """

    def __init__(self, name, genome_dir=None):
        self.genome_dir = get_genome_dir(genome_dir)
        self.name, self.filename = self._get_name_and_filename(name)
        super(Genome, self).__init__(self.filename)

        metadata = self._read_metadata()
        self.tax_id = metadata.get("tax_id")
        self.assembly_accession = metadata.get("assembly_accession")
        self._gap_sizes = None
        self.props = {}

        for plugin in get_active_plugins():
            self.props[plugin.name()] = plugin.get_properties(self)

    def _get_name_and_filename(self, name):
        """
        name can be just a name (e.g. hg38) or an abspath to a fasta file or the fasta's folder

        returns the name and the abspath to the (closest) fasta file
        """
        stripped_name = get_localname(name)
        stripped_name = re.sub(".fa(.gz)?$", "", stripped_name)

        if os.path.isfile(name):
            filename = name
            name = os.path.basename(stripped_name)
        elif os.path.isdir(name) and glob_ext_files(name)[0].startswith(
            os.path.basename(name) + ".fa"
        ):
            filename = glob_ext_files(name)[0]
            name = os.path.basename(stripped_name)
        else:
            fasta_dir = os.path.join(self.genome_dir, stripped_name)
            filenames = glob_ext_files(fasta_dir)
            if len(filenames) == 1:
                filename = filenames[0]
            elif len(filenames) == 0:
                raise FileNotFoundError(
                    f"no *.fa files found in genome_dir {fasta_dir}"
                )
            else:
                filename = os.path.join(fasta_dir, stripped_name + ".fa")
                if filename not in filenames:
                    filename += ".gz"
                if filename not in filenames:
                    raise Exception(
                        f"Multiple fasta files found, but not {stripped_name}.fa!"
                    )
            name = stripped_name

        return [name, filename]

    def _read_metadata(self):
        """
        Read genome metadata from genome README.txt (if it exists).
        """
        readme = os.path.join(self.genome_dir, self.name, "README.txt")
        metadata, lines = read_readme(readme)

        update_metadata = False
        if metadata["genome url"] != "na":
            for key in ["provider", "tax_id", "assembly_accession"]:
                if key not in metadata:
                    update_metadata = True
                    break
        if not update_metadata:
            return metadata

        print(f"Updating metadata in README.txt", file=sys.stderr)

        if "provider" not in metadata:
            url = metadata.get("genome url", "").lower()
            for provider in ["Ensembl", "UCSC", "NCBI"]:
                if provider.lower() in url:
                    metadata["provider"] = provider
            else:
                metadata["provider"] = "Unknown"

        if "tax_id" not in metadata or "assembly_accession" not in metadata:
            name = safe(metadata.get("original name", ""))
            p = None
            genome = None
            if metadata["provider"] != "Unknown":
                p = ProviderBase.create(metadata["provider"])
                if name and name in p.genomes:
                    genome = p.genomes[name]

            if "tax_id" not in metadata:
                metadata["tax_id"] = "na"
                if genome:
                    taxid = p.genome_taxid(genome)
                    taxid = str(taxid) if taxid != 0 else "na"
                    metadata["tax_id"] = taxid

            if "assembly_accession" not in metadata:
                metadata["assembly_accession"] = "na"
                if genome:
                    accession = p.assembly_accession(genome)
                    metadata["assembly_accession"] = accession

        write_readme(readme, metadata, lines)
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

    def gap_sizes(self):
        """Return gap sizes per chromosome.

        Returns
        -------
        gap_sizes : dict
            a dictionary with chromosomes as key and the total number of
            Ns as values
        """
        if not self._gap_sizes:
            gap_file = self.filename + ".sizes"

            # generate gap file if not found
            if not os.path.exists(gap_file):
                generate_gap_bed(self.filename, gap_file)

            self._gap_sizes = {}
            with open(gap_file) as f:
                for line in f:
                    chrom, start, end = line.strip().split("\t")
                    start, end = int(start), int(end)
                    self._gap_sizes[chrom] = self._gap_sizes.get(chrom, 0) + end - start
        return self._gap_sizes

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
        if self._gap_sizes is None:
            self.gap_sizes()

        # dict of chromosome sizes after subtracting the number of Ns
        sizes = dict(
            [
                (chrom, len(self[chrom]) - self._gap_sizes.get(chrom, 0))
                for chrom in chroms
            ]
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
