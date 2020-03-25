import os.path
import re
import sys

from collections.abc import Iterable
from bisect import bisect
from pyfaidx import Fasta, Sequence
from random import random

from genomepy.exceptions import GenomeDownloadError
from genomepy.plugin import get_active_plugins
from genomepy.provider import ProviderBase
from genomepy.utils import (
    get_genome_dir,
    get_localname,
    glob_ext_files,
    generate_gap_bed,
)


def _weighted_selection(l, n):
    """
    Selects  n random elements from a list of (weight, item) tuples.
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
        self.name, filename = self._get_name_and_filename(name)
        super(Genome, self).__init__(filename)

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
        stripped_name = re.sub(".gz$", "", stripped_name)
        stripped_name = re.sub(".fa$", "", stripped_name)

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
        metadata = {}
        readme = os.path.join(self.genome_dir, self.name, "README.txt")
        if os.path.exists(readme):
            with open(readme) as f:
                metadata = {}
                for line in f.readlines():
                    vals = line.strip().split(":")
                    metadata[vals[0].strip()] = (":".join(vals[1:])).strip()

            update_metadata = False
            for key in ["provider", "tax_id", "assembly_accession"]:
                if key not in metadata:
                    update_metadata = True
                    break
            if not update_metadata:
                return metadata

            print(f"Updating metadata in README.txt", file=sys.stderr)

            if "provider" not in metadata:
                if "ensembl" in metadata.get("url", ""):
                    metadata["provider"] = "Ensembl"
                elif "ucsc" in metadata.get("url", ""):
                    metadata["provider"] = "UCSC"
                elif "ncbi" in metadata.get("url", ""):
                    metadata["provider"] = "NCBI"
                else:
                    metadata["provider"] = "URL"

            if metadata.get("provider", "").lower() in ["ensembl", "ucsc", "ncbi"]:
                if "tax_id" not in metadata or "assembly_accession" not in metadata:
                    p = ProviderBase.create(metadata["provider"])

                    if "tax_id" not in metadata:
                        try:
                            metadata["tax_id"] = p.genome_taxid(
                                metadata["original name"]
                            )
                        except GenomeDownloadError:
                            print(
                                f"Could not update tax_id of {self.name}",
                                file=sys.stderr,
                            )
                    if "assembly_accession" not in metadata:
                        try:
                            metadata["assembly_accession"] = p.assembly_accession(
                                metadata["original name"]
                            )
                        except GenomeDownloadError:
                            print(
                                f"Could not update assembly_accession of {self.name}",
                                file=sys.stderr,
                            )

            with open(readme, "w") as f:
                for k, v in metadata.items():
                    print(f"{k}: {v}", file=f)
        return metadata

    def _bed_to_seqs(self, track, stranded=False, extend_up=0, extend_down=0):
        bufsize = 10000
        with open(track) as fin:
            lines = fin.readlines(bufsize)
            while lines:
                for line in lines:
                    if line.startswith("#") or line.startswith("track"):
                        continue

                    vals = line.strip().split("\t")
                    try:
                        start, end = int(vals[1]), int(vals[2])
                    except ValueError:
                        raise

                    rc = False
                    if stranded:
                        try:
                            rc = vals[5] == "-"
                        except IndexError:
                            pass

                    starts = [start]
                    ends = [end]

                    chrom = vals[0]

                    # BED12
                    if len(vals) == 12:
                        starts = [int(x) for x in vals[11].split(",")[:-1]]
                        sizes = [int(x) for x in vals[10].split(",")[:-1]]
                        starts = [start + x for x in starts]
                        ends = [start + size for start, size in zip(starts, sizes)]
                    name = "{}:{}-{}".format(chrom, start, end)
                    try:
                        name = " ".join((name, vals[3]))
                    except Exception:
                        pass

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

                lines = fin.readlines(bufsize)

    def _region_to_seqs(self, track, extend_up=0, extend_down=0):
        bufsize = 10000
        if isinstance(track, list):
            for name in track:
                chrom, coords = name.split(":")
                start, end = [int(c) for c in coords.split("-")]
                start += 1
                start -= extend_up
                end += extend_down
                seq = self.get_seq(chrom, start, end)
                yield Sequence(name, seq.seq)
        else:
            with open(track) as fin:
                lines = fin.readlines(bufsize)
                while lines:
                    for line in lines:
                        name = line.strip()
                        chrom, coords = name.split(":")
                        start, end = [int(c) for c in coords.split("-")]
                        start += 1
                        start -= extend_up
                        end += extend_down
                        seq = self.get_seq(chrom, start, end)
                        yield Sequence(name, seq.seq)

                    lines = fin.readlines(bufsize)

    @staticmethod
    def get_track_type(track):
        region_p = re.compile(r"^(.+):(\d+)-(\d+)$")
        if not isinstance(track, (str, bytes)) and isinstance(track, Iterable):
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
        track_type = self.get_track_type(track)
        if track_type == "interval":
            seqqer = self._region_to_seqs(
                track, extend_up=extend_up, extend_down=extend_down
            )
        else:
            seqqer = self._bed_to_seqs(
                track, stranded=stranded, extend_up=extend_up, extend_down=extend_down
            )

        if fastafile:
            with open(fastafile, "w") as fout:
                for seq in seqqer:
                    fout.write("{}\n".format(seq.__repr__()))
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
            # gap_file = self.props["gaps"]["gaps"]
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

    def get_random_sequences(self, n=10, length=200, chroms=None, max_n=0.1):
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

        Returns
        -------
        coords : list
            List with [chrom, start, end] genomic coordinates.
        """
        retries = 100
        cutoff = length * max_n
        if not chroms:
            chroms = self.keys()

        if self._gap_sizes is None:
            self.gap_sizes()
        sizes = dict(
            [
                (chrom, len(self[chrom]) - self._gap_sizes.get(chrom, 0))
                for chrom in chroms
            ]
        )

        lengths = [
            (sizes[x], x)
            for x in chroms
            if sizes[x] / len(self[x]) > 0.1 and sizes[x] > 10 * length
        ]
        chroms = _weighted_selection(lengths, n)
        coords = []

        count = {}
        for chrom in chroms:
            if chrom in count:
                count[chrom] += 1
            else:
                count[chrom] = 1

        for chrom in chroms:
            for _ in range(retries):
                start = int(random() * (sizes[chrom] - length))
                end = start + length
                count_n = self[chrom][start:end].seq.upper().count("N")
                if count_n <= cutoff:
                    break
            else:
                raise Exception("Unexpected error")
            if count_n > cutoff:
                raise ValueError(
                    "Failed to find suitable non-N sequence for {}".format(chrom)
                )

            coords.append([chrom, start, end])

        return coords
