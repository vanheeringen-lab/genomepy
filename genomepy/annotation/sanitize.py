from typing import Optional

import pandas as pd
from loguru import logger

from genomepy.annotation.utils import _check_property, bed_from_gtf, write_annot
from genomepy.files import _open, update_readme


def sanitize(
    self,
    match_contigs: Optional[bool] = True,
    filter_contigs: Optional[bool] = True,
):
    """
    Compares the genome and gene annotations.
    If the contig names do not conform, attempt to fix this.
    If the contig names conform (after fixing), filter the annotation contigs for genome contigs.
    Log the results in the README.

    Parameters
    ----------
    self: Annotation class instance

    match_contigs: bool, optional
        attempt to fix contig names?

    filter_contigs: bool, optional
        remove contigs from the annotations that are missing from the genome?
    """
    if self.genome_file is None:
        logger.error("A genome is required for sanitizing!")
        return
    _check_property(self.readme_file, "README.txt")

    extra_lines = []
    if is_overlapping(self):
        status = "contigs match but not filtered"
        if filter_contigs:
            status = "contigs match"
            contigs_filtered_out = filter_genome_contigs(self)
            if contigs_filtered_out:
                bed_from_gtf(self.annottation_gtf_file, self.annotation_bed_file)
                status = "contigs match and filtered"
                extra_lines = [
                    "",
                    "The following contigs were filtered out of the gene annotation:",
                    f"{', '.join(contigs_filtered_out)}",
                ]
        update_readme(self.readme_file, {"sanitized annotation": status}, extra_lines)
        return

    if match_contigs is False:
        return

    matching_contig_index = conforming_index(self)
    if is_conformable(matching_contig_index) is False:
        # example not possible: ASM2732v1
        update_readme(self.readme_file, {"sanitized annotation": "not possible"})
        return

    conversion_dict, duplicate_contigs = contig_conversion_dict(
        self.genome_file, matching_contig_index
    )
    if duplicate_contigs:
        logger.warning(
            "\nThe genome contains duplicate contig names, this should really not happen!\n"
            "You may wish to consider filtering these out, or look for another genome assembly..."
            "The following duplicate contigs were found: "
            f"{', '.join(duplicate_contigs)}.\n"
        )
    missing_contigs = conform_gtf(
        self.annotation_gtf_file, conversion_dict, filter_contigs
    )
    bed_from_gtf(self.annotation_gtf_file, self.annotation_bed_file)
    status = "contigs fixed"
    if missing_contigs:
        status = "contigs fixed and filtered"
        extra_lines = [
            "",
            "The following contigs were filtered out of the gene annotation:",
            f"{', '.join(missing_contigs)}",
        ]
        if filter_contigs is False:
            status = "contigs fixed but not filtered"
            extra_lines[1] = (
                "The following contigs could not be sanitized"
                " in the gene annotation, and were kept as-is:"
            )
        logger.info(" ".join(extra_lines))
    update_readme(self.readme_file, {"sanitized annotation": status}, extra_lines)


def is_overlapping(self):
    """
    Check if genome and annotation contigs conform.
    Returns bool.
    """
    overlapping_contigs = set(self.genome_contigs) & set(self.annotation_contigs)
    return bool(overlapping_contigs)


def filter_genome_contigs(self):
    """
    filter the gtf for contigs in the genome.
    return a list of discarded contigs.
    """
    old = pd.read_csv(self.annotation_gtf_file, sep="\t", header=None, comment="#")
    filter_func = old[0].isin(self.genome_contigs)
    new = old[filter_func]
    write_annot(new, self.annotation_gtf_file)

    discarded_contigs = list(set(old[~filter_func][0]))
    return discarded_contigs


def conforming_index(self) -> int:
    """
    Check if the annotation contigs can be made to conform the genome contigs.
    Returns int: -1 = not conformable, 0 = conforming, >0 = index of conforming element.
    """
    # example conformable files
    # genome.fa:
    #
    # >chr1 alternate_name_chr1 more_info                    <-- header
    # ATCGATCGATCGATCG                                       <-- sequence
    #
    # genome.annotation.gtf:
    #
    # seqname              source    feature start  end ...  <-- header
    # alternate_name_chr1  provider  exon    0      100 ...
    #
    # returns 1 because genome header[1] matches the GTF contig names

    # all headers in the file should be formatted the same. so take the first.
    header = []
    with _open(self.genome_file, "r") as fa:
        for line in fa:
            if line.startswith(">"):
                header = line.strip(">\n").split(" ")
                break

    for contig in self.annotation_contigs:
        if contig in header:
            matching_contig_index = header.index(contig)
            return matching_contig_index
    return -1


def is_conformable(matching_contig_index: int):
    if matching_contig_index == 0:
        raise ValueError("genome and annotation are conforming.")
    return matching_contig_index > 0


def contig_conversion_dict(genome_file, matching_element: int):
    """
    build a conversion dict
    returns the dict and a list of duplicate contig names (should be empty!)
    """
    conversion_table = {}
    duplicate_contigs = []
    with _open(genome_file, "r") as fa:
        for line in fa:
            if line.startswith(">"):
                line = line.strip(">\n").split(" ")
                if line[matching_element] in conversion_table:
                    duplicate_contigs.append(line[matching_element])
                conversion_table[line[matching_element]] = line[0]
    return conversion_table, list(set(duplicate_contigs))


def conform_gtf(gtf_file, conversion_dict: dict, filter_contigs: Optional[bool] = True):
    """
    Create an annotation gtf file with contigs conforming to the genome
    Overwrites the existing gtf file.
    Returns a list of contigs not present in the genome
    """
    old = pd.read_csv(gtf_file, sep="\t", header=None, comment="#")
    converted_contigs = old[0].map(conversion_dict)  # contigs missing become NaNs
    nans = converted_contigs.isna()

    new = old.copy()
    new[0] = converted_contigs.fillna(old[0])  # convert contigs (keeps missing)
    if filter_contigs:
        new = new[~nans]  # remove missing contigs
    write_annot(new, gtf_file)

    missing_contigs = list(set(old[nans][0]))
    return missing_contigs
