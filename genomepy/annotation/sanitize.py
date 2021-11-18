"""Annotation.sanitize method"""
from itertools import compress

from genomepy.annotation.utils import _check_property, write_annot
from genomepy.files import _open, update_readme


def _sanitize(self, match=True, filter=True, overwrite=False):  # noqa
    """
    Match the contigs names of the gene annotations to the genome's.

    First, match the contig names if possible.
    Second, remove contig names not found in the genome.
    Third, save the results and log this in the README.

    Parameters
    ----------
    match: bool, optional
        match annotation contig names to the genome contig names (default is True)
    filter: bool, optional
        remove annotation contig names not found in the genome contig names (default is True)
    overwrite: bool, optional
        update the annotation files on disk, and log this in the README (default is False).

    Returns
    -------
    Annotation class
        updated attributes
    """
    _check_property(self.genome_file, f"{self.name}.fa")

    cd = {}
    if match:
        cd = _match_contigs(self)

    mc = []
    if filter:
        mc = _filter_contigs(self)

    if overwrite:
        write_annot(self.gtf, self.annotation_gtf_file)
        write_annot(self.bed, self.annotation_bed_file)
        _document_sanitizing(self, cd, mc)


def _match_contigs(self):
    """
    Sometimes, the gene annotation contig names are different from the genome contig names.
    If the genome contains multiple headers, these may contain the annotation contig names.

    Match the gtf contig names to the genome contig names where possible.
    """
    # are there contigs in the annotation not found in the genome?
    missing_contigs = bool(set(self.annotation_contigs).difference(self.genome_contigs))

    if missing_contigs is False:
        return  # nothing to do

    # get the full headers from the fasta (e.g. '>chr1 chromosome1 note1 note2')
    headers = _full_genome_headers(self.genome_file)

    # filter for multi-head headers (whitespace delimited)
    multi_headed_contigs = [" " in h or "\t" in h for h in headers]

    if not any(multi_headed_contigs):
        return  # nothing to work with

    headers = [h.split() for h in compress(headers, multi_headed_contigs)]
    conversion_dict = _contig_conversion_dict(headers, self.annotation_contigs)

    self.gtf = self.gtf.replace({"seqname": conversion_dict})
    self.bed = self.bed.replace({"chrom": conversion_dict})
    return conversion_dict


def _filter_contigs(self):
    """
    Some tools throw a fit when the gene annotation contains
    contigs that are missing from the genome.

    Remove contigs not found in the genome.
    """
    # set of contigs in the genome that are not in the annotation
    missing_contigs = set(self.annotation_contigs).difference(self.genome_contigs)

    self.gtf = self.gtf[~self.gtf.seqname.isin(missing_contigs)]
    self.bed = self.bed[~self.bed.chrom.isin(missing_contigs)]
    return missing_contigs


def _document_sanitizing(self, cd, mc):
    _check_property(self.readme_file, "README.txt")
    status = ""
    extra_lines = []

    matched = bool(cd)
    if matched:
        status += f"{len(cd)} contigs were renamed. "

    filtered = bool(mc)
    if filtered:
        status += f"{len(mc)} contigs were removed (see below)."
        extra_lines = [
            "",
            "The following contigs were filtered out of the gene annotation:",
            f"{', '.join(mc)}",
        ]

    if not matched and not filtered:
        status = "no changes made."

    update_readme(self.readme_file, {"sanitized annotation": status}, extra_lines)


def _full_genome_headers(genome_file):
    headers = []
    with _open(genome_file) as fa:
        for line in fa:
            if line[0] == ">":
                headers.append(line.strip("> \n"))
    return headers


def _contig_conversion_dict(genome_headers, annotation_contigs):
    conversion_table = {}
    for header in genome_headers:
        for contig in annotation_contigs:
            if contig in header:
                conversion_table[contig] = header[0]
                break  # next header
    return conversion_table
