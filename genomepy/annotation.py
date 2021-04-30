import csv
import os
import re
import subprocess as sp
import sys
from tempfile import mkdtemp
from typing import Optional, Union

import pandas as pd

from genomepy.utils import (
    get_genomes_dir,
    gzip_and_name,
    gunzip_and_name,
    glob_ext_files,
    update_readme,
    rm_rf,
    _open,
)


class Annotation:
    """
    Gene annotation related functions

    Parameters
    ----------
    name : str
        Genome name

    genomes_dir : str
        Genomes installation directory

    Returns
    -------
    Annotation object
    """

    def __init__(self, name: str, genomes_dir: str = None):
        self.name = name
        self.genome_dir = os.path.join(get_genomes_dir(genomes_dir), name)
        self.readme_file = os.path.join(self.genome_dir, "README.txt")
        self.annotation_gtf_file = self._get_genome_file("annotation.gtf")
        self.annotation_bed_file = self._get_genome_file("annotation.bed")

        # variables that may be None
        self.genome_file = self._get_genome_file("fa", check_exists=False)
        self.index_file = self.genome_file + ".fai" if self.genome_file else None
        self.sizes_file = self.genome_file + ".sizes" if self.genome_file else None

        # variables to store if called
        self.genome_contigs = []
        self.annotation_contigs = []

    def _get_genome_file(self, ext: str, check_exists: Optional[bool] = True):
        """
        Returns the filepath to a single (gzipped) file in the genome_dir with matching ext.
        """
        ext_files = glob_ext_files(self.genome_dir, ext)
        fname = [f for f in ext_files if f"{self.name}.{ext}" in f]
        if len(fname) == 1:
            filepath = os.path.join(self.genome_dir, fname[0])
            return filepath

        # error handling
        if len(fname) > 1:
            raise IndexError(
                f"Too many ({len(fname)}) files match '{self.name}.{ext}(.gz)' "
                f"in directory {self.genome_dir}"
            )
        if check_exists:
            raise FileNotFoundError(
                f"could not find '{self.name}.{ext}(.gz)' in directory {self.genome_dir}"
            )
        return None

    @property
    def genome_contigs(self):
        return self.__genome_contigs

    @genome_contigs.setter
    def genome_contigs(self, genome_contig_list):
        self.__genome_contigs = genome_contig_list

    @genome_contigs.getter
    def genome_contigs(self):
        if not self.__genome_contigs:
            return get_column(self.sizes_file)

    @property
    def annotation_contigs(self):
        return self.__annotation_contigs

    @annotation_contigs.setter
    def annotation_contigs(self, annotation_contig_list):
        self.__annotation_contigs = annotation_contig_list

    @annotation_contigs.getter
    def annotation_contigs(self):
        if not self.__annotation_contigs:
            return list(set(get_column(self.annotation_bed_file)))

    @staticmethod
    def filter_regex(
        annotation_file: str,
        regex: Optional[str] = ".*",
        invert_match: Optional[bool] = False,
    ):
        """
        Filter annotation file (works with both gtf and bed) by contig name.

        annotation_file: path to annotation file.
        regex: regex string to keep.
        invert_match: keep contigs NOT matching the regex string.

        return a list of discarded contigs.
        """
        old = pd.read_csv(annotation_file, sep="\t", header=None, comment="#")
        pattern = re.compile(regex)
        filter_func = old[0].map(lambda x: bool(pattern.match(x)) is not invert_match)
        new = old[filter_func]
        new.to_csv(
            annotation_file, sep="\t", header=None, index=None, quoting=csv.QUOTE_NONE
        )

        discarded_contigs = list(set(old[~filter_func][0]))
        return discarded_contigs

    def filter_genome_contigs(self):
        """
        filter the gtf for contigs in the genome.
        return a list of discarded contigs.
        """
        old = pd.read_csv(self.annotation_gtf_file, sep="\t", header=None, comment="#")
        filter_func = old[0].isin(self.genome_contigs)
        new = old[filter_func]
        new.to_csv(
            self.annotation_gtf_file,
            sep="\t",
            header=None,
            index=None,
            quoting=csv.QUOTE_NONE,
        )

        discarded_contigs = list(set(old[~filter_func][0]))
        return discarded_contigs

    def _is_conforming(self):
        """
        Check if genome and annotation contigs conform.
        Returns bool.
        """
        overlapping_contigs = set(self.genome_contigs) & set(self.annotation_contigs)
        return bool(overlapping_contigs)

    def _conforming_index(self) -> int:
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

    @staticmethod
    def _is_conformable(matching_contig_index: int):
        if matching_contig_index == 0:
            raise ValueError("genome and annotation are conforming.")
        return matching_contig_index > 0

    def _contig_conversion_dict(self, matching_element: int):
        """
        build a conversion dict
        returns the dict and a list of duplicate contig names (should be empty!)
        """
        conversion_table = {}
        duplicate_contigs = []
        with _open(self.genome_file, "r") as fa:
            for line in fa:
                if line.startswith(">"):
                    line = line.strip(">\n").split(" ")
                    if line[matching_element] in conversion_table:
                        duplicate_contigs.append(line[matching_element])
                    conversion_table[line[matching_element]] = line[0]
        return conversion_table, list(set(duplicate_contigs))

    def _conform_gtf(
        self, conversion_dict: dict, filter_contigs: Optional[bool] = True
    ):
        """
        Create an annotation gtf file with contigs conforming to the genome
        Overwrites the existing gtf file.
        Returns a list of contigs not present in the genome
        """
        old = pd.read_csv(self.annotation_gtf_file, sep="\t", header=None, comment="#")
        converted_contigs = old[0].map(conversion_dict)  # contigs missing become NaNs
        nans = converted_contigs.isna()

        new = old.copy()
        new[0] = converted_contigs.fillna(old[0])  # convert contigs (keeps missing)
        if filter_contigs:
            new = new[~nans]  # remove missing contigs
        new.to_csv(
            self.annotation_gtf_file,
            sep="\t",
            header=None,
            index=None,
            quoting=csv.QUOTE_NONE,
        )

        missing_contigs = list(set(old[nans][0]))
        return missing_contigs

    def bed_from_gtf(self):
        """
        Create an annotation bed file from an annotation gtf file.
        Overwrites the existing bed file.
        """
        tmp_dir = mkdtemp(dir=self.genome_dir)

        # unzip if needed
        self.annotation_gtf_file, is_unzipped = gunzip_and_name(
            self.annotation_gtf_file
        )
        new_bed_file = os.path.join(
            tmp_dir, os.path.basename(re.sub(r"\.gz$", "", self.annotation_bed_file))
        )

        cmd = "gtfToGenePred {0} /dev/stdout | genePredToBed /dev/stdin {1}"
        sp.check_call(cmd.format(self.annotation_gtf_file, new_bed_file), shell=True)

        # gzip if needed
        self.annotation_gtf_file = gzip_and_name(self.annotation_gtf_file, is_unzipped)
        new_bed_file = gzip_and_name(
            new_bed_file, self.annotation_bed_file.endswith(".gz")
        )

        os.replace(new_bed_file, self.annotation_bed_file)
        rm_rf(tmp_dir)

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

        match_contigs: attempt to fix contig names?
        filter_contigs: remove contigs from the annotations that are missing from the genome?
        """
        if self.genome_file is None:
            sys.stderr.write("A genome is required for sanitizing!")
            return

        extra_lines = []
        if self._is_conforming():
            status = "contigs match but not filtered"
            if filter_contigs:
                status = "contigs match"
                contigs_filtered_out = self.filter_genome_contigs()
                if contigs_filtered_out:
                    self.bed_from_gtf()
                    status = "contigs match and filtered"
                    extra_lines = [
                        "",
                        "The following contigs were filtered out of the gene annotation:",
                        f"{', '.join(contigs_filtered_out)}",
                    ]
            update_readme(
                self.readme_file, {"sanitized annotation": status}, extra_lines
            )
            return

        if match_contigs is False:
            return

        matching_contig_index = self._conforming_index()
        if self._is_conformable(matching_contig_index) is False:
            # example not possible: ASM2732v1
            update_readme(self.readme_file, {"sanitized annotation": "not possible"})
            return

        conversion_dict, duplicate_contigs = self._contig_conversion_dict(
            matching_contig_index
        )
        if duplicate_contigs:
            sys.stderr.write(
                "\nThe genome contains duplicate contig names, this should really not happen!\n"
                "You may wish to consider filtering these out, or look for another genome assembly..."
                "The following duplicate contigs were found: "
                f"{', '.join(duplicate_contigs)}.\n"
            )
        missing_contigs = self._conform_gtf(conversion_dict, filter_contigs)
        self.bed_from_gtf()
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
            sys.stderr.write(" ".join(extra_lines))
        update_readme(self.readme_file, {"sanitized annotation": status}, extra_lines)


def get_column(
    fname: str,
    n: Optional[int] = 0,
    sep: str = "\t",
    comment_char: Union[None, str] = "#",
):
    """
    Return the nth column from a separated file.
    """

    def comment(_: str):
        return False

    if isinstance(comment_char, str):

        def comment(string: str):  # noqa: F811
            return string.startswith(comment_char)

    column = []
    with _open(fname) as f:
        for line in f:
            line = str(line).strip()
            if not comment(line):
                column.append(line.split(sep)[n])
    return column
