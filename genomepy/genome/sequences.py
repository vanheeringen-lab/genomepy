import re
from bisect import bisect
from random import random

from pyfaidx import Sequence


def track2fasta(
    self, track, fastafile=None, stranded=False, extend_up=0, extend_down=0
):
    """
    Return a list of fasta sequences as Sequence objects
    as directed from the track(s).

    Parameters
    ----------
    self: Genome class instance

    track: list/region file/bed file
        region(s) you wish to translate to fasta.
        Example input files can be found in genomepy/tests/data/regions.*

    fastafile: bool , optional
        return Sequences as list or save to file? (default: list)

    stranded: bool , optional
        return sequences for both strands? Required BED6 (or higher) as input (default: False)

    extend_up: int , optional
        extend the sequences up? (command is strand sensitive, default: 0)

    extend_down: int , optional
        extend the sequences down? (command is strand sensitive, default: 0)
    """
    track_type = get_track_type(track)
    if track_type == "interval":
        seqqer = regions_to_seqs(
            self, track, extend_up=extend_up, extend_down=extend_down
        )
    else:
        seqqer = bed_to_seqs(
            self, track, stranded=stranded, extend_up=extend_up, extend_down=extend_down
        )

    if fastafile:
        with open(fastafile, "w") as fout:
            for seq in seqqer:
                fout.write(f"{seq.__repr__()}\n")
    else:
        return [seq for seq in seqqer]


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


def regions_to_seqs(self, track, extend_up=0, extend_down=0):
    if isinstance(track, list):
        for region in track:
            name = region.strip()
            seq = region_to_seq(self, name, extend_up, extend_down)
            yield Sequence(name, seq)
    else:
        with open(track) as fin:
            bufsize = 10000
            lines = fin.readlines(bufsize)
            for region in lines:
                name = region.strip()
                seq = region_to_seq(self, name, extend_up, extend_down)
                yield Sequence(name, seq)

                # load more lines if needed
                lines += fin.readlines()


def region_to_seq(self, region, extend_up=0, extend_down=0):
    chrom, coords = region.strip().split(":")
    start, end = [int(c) for c in coords.split("-")]
    start += 1
    start -= extend_up
    end += extend_down
    seq = self.get_seq(chrom, start, end)
    return seq.seq


def bed_to_seqs(self, track, stranded=False, extend_up=0, extend_down=0):
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


def get_random_sequences(
    self, n=10, length=200, chroms=None, max_n=0.1, outtype="list"
):
    """
    Return random genomic sequences.

    Parameters
    ----------
    self: Genome class instance

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
    chroms = weighted_selection(lengths, n)

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


def weighted_selection(list_of_tuples, n):
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
