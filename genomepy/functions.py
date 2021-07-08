"""Module-level functions."""
import os
import re
from functools import singledispatch
from io import TextIOWrapper
from tempfile import NamedTemporaryFile
from typing import Optional

import pyfaidx
from appdirs import user_config_dir
from Bio.SeqIO.FastaIO import SimpleFastaParser

from genomepy.annotation import Annotation
from genomepy.config import config
from genomepy.exceptions import GenomeDownloadError
from genomepy.files import (
    _apply_fasta_regex_func,
    bgzip_and_name,
    glob_ext_files,
    gzip_and_name,
    read_readme,
    update_readme,
)
from genomepy.genome import Genome
from genomepy.online import check_url
from genomepy.plugins import get_active_plugins
from genomepy.providers import download_assembly_report, online_providers
from genomepy.utils import (
    get_genomes_dir,
    get_localname,
    mkdir_p,
    rm_rf,
    safe,
    try_except_pass,
)


def list_available_genomes(provider=None):
    """
    List all available genomes.

    Parameters
    ----------
    provider : str, optional
        List genomes from specific provider. Genomes from all
        providers will be returned if not specified.

    Returns
    -------
    list with genome names
    """
    for p in online_providers(provider):
        for row in p.list_available_genomes():
            yield list(row[:1]) + [p.name] + list(row[1:])


def list_installed_genomes(genomes_dir: str = None):
    """
    List all locally available genomes.

    Parameters
    ----------
    genomes_dir : str
        Directory with genomes installed by genomepy.

    Returns
    -------
    list with genome names
    """
    genomes_dir = get_genomes_dir(genomes_dir, check_exist=False)
    if os.path.exists(genomes_dir):
        return [
            subdir
            for subdir in os.listdir(genomes_dir)
            if _is_genome_dir(os.path.join(genomes_dir, subdir))
        ]
    return []


def install_genome(
    name: str,
    provider: Optional[str] = None,
    genomes_dir: Optional[str] = None,
    localname: Optional[str] = None,
    mask: Optional[str] = "soft",
    keep_alt: Optional[bool] = False,
    regex: Optional[str] = None,
    invert_match: Optional[bool] = False,
    bgzip: Optional[bool] = None,  # None -> check config. False -> dont check.
    annotation: Optional[bool] = False,
    only_annotation: Optional[bool] = False,
    skip_matching: Optional[bool] = False,
    skip_filter: Optional[bool] = False,
    threads: Optional[int] = 1,
    force: Optional[bool] = False,
    **kwargs: Optional[dict],
):
    """
    Install a genome.

    Parameters
    ----------
    name : str
        Genome name

    provider : str , optional
        Provider name. will try Ensembl, UCSC and NCBI (in that order) if not specified.

    genomes_dir : str , optional
        Where to store the fasta files

    localname : str , optional
        Custom name for this genome.

    mask : str , optional
        Default is 'soft', choices 'hard'/'soft/'none' for respective masking level.

    keep_alt : bool , optional
        Some genomes contain alternative regions. These regions cause issues with
        sequence alignment, as they are inherently duplications of the consensus regions.
        Set to true to keep these alternative regions.

    regex : str , optional
        Regular expression to select specific chromosome / scaffold names.

    invert_match : bool , optional
        Set to True to select all chromosomes that don't match the regex.

    bgzip : bool , optional
        If set to True the genome FASTA file will be compressed using bgzip,
        and gene annotation will be compressed with gzip.
        If not specified, the setting from the configuration file will be used.

    threads : int , optional
        Build genome index using multithreading (if supported). Default: lowest of 8/all threads

    force : bool , optional
        Set to True to overwrite existing files.

    annotation : bool , optional
        If set to True, download gene annotation in BED and GTF format.

    only_annotation : bool , optional
        If set to True, only download the annotation files.

    skip_matching : bool , optional
        If set to True, contigs in the annotation not matching
        those in the genome will not be corrected.

    skip_filter : bool , optional
        If set to True, the gene annotations will not be filtered to match the genome contigs.

    kwargs : dict , optional
        Provider specific options.
        toplevel : bool , optional
            Ensembl only: Always download the toplevel genome. Ignores potential primary assembly.

        version : int , optional
            Ensembl only: Specify release version. Default is latest.

        to_annotation : text , optional
            URL only: direct link to annotation file.
            Required if this is not the same directory as the fasta.
    """
    name = safe(name)
    localname = get_localname(name, localname)
    genomes_dir = get_genomes_dir(genomes_dir, check_exist=False)
    out_dir = os.path.join(genomes_dir, localname)
    genome_file = os.path.join(out_dir, f"{localname}.fa")
    provider = _provider_selection(name, localname, genomes_dir, provider)

    # check which files need to be downloaded
    genome_found = _is_genome_dir(out_dir)
    download_genome = (
        genome_found is False or force is True
    ) and only_annotation is False
    annotation_found = bool(glob_ext_files(out_dir, "annotation.gtf")) and bool(
        glob_ext_files(out_dir, "annotation.bed")
    )
    download_annotation = (annotation_found is False or force is True) and any(
        [
            annotation,
            only_annotation,
            skip_matching,
            skip_filter,
            kwargs.get("to_annotation"),
            kwargs.get("ucsc_annotation_type"),
        ]
    )

    genome = None
    genome_downloaded = False
    if download_genome:
        if force:
            _delete_extensions(out_dir, ["fa", "fai"])
        provider.download_genome(
            name,
            genomes_dir,
            mask=mask,
            localname=localname,
            **kwargs,
        )
        genome_found = True
        genome_downloaded = True

        # Filter genome
        _filter_genome(genome_file, regex, invert_match, keep_alt)

        # Generates a Fasta object and the genome index, gaps and sizes files
        genome = Genome(localname, genomes_dir=genomes_dir)

        # Download the NCBI assembly report
        asm_report = os.path.join(out_dir, "assembly_report.txt")
        asm_acc = genome.assembly_accession
        if not os.path.exists(asm_report) and asm_acc != "na":
            download_assembly_report(asm_acc, asm_report)

        # Export installed genome(s)
        generate_env(genomes_dir=genomes_dir)

    annotation_downloaded = False
    if download_annotation:
        if force:
            _delete_extensions(out_dir, ["annotation.gtf", "annotation.bed"])
        provider.download_annotation(name, genomes_dir, localname=localname, **kwargs)
        annotation_downloaded = bool(
            glob_ext_files(out_dir, "annotation.gtf")
        ) and bool(glob_ext_files(out_dir, "annotation.bed"))

    if annotation_downloaded:
        annotation = Annotation(localname, genomes_dir)
        if genome_found and not (skip_matching and skip_filter):
            annotation.sanitize(not skip_matching, not skip_filter, True)

    # Run active plugins (also if the genome was downloaded earlier)
    if genome_found:
        genome = genome if genome else Genome(localname, genomes_dir=genomes_dir)
        for plugin in get_active_plugins():
            plugin.after_genome_download(genome, threads, force)

    # zip files downloaded now
    if bgzip is True or (bgzip is None and config.get("bgzip")):
        if genome_downloaded:
            bgzip_and_name(genome.filename)
        if annotation_downloaded:
            gzip_and_name(annotation.annotation_gtf_file)
            gzip_and_name(annotation.annotation_bed_file)

    return genome


def generate_env(fname: str = "exports.txt", genomes_dir: str = None):
    """
    Generate file with exports.

    By default this is .config/genomepy/exports.txt.

    An alternative file name or file path is accepted too.

    Parameters
    ----------
    fname: str, optional
        Absolute path or name of the output file.

    genomes_dir: str, optional
        Directory with installed genomes to export.
    """
    fname1 = os.path.expanduser(fname)
    fname2 = os.path.join(user_config_dir("genomepy"), fname)
    fname = fname1 if os.path.isabs(fname1) else fname2
    mkdir_p(os.path.dirname(fname))
    with open(fname, "w") as fout:
        for env in _generate_exports(genomes_dir):
            fout.write(f"{env}\n")


def _generate_exports(genomes_dir: str = None):
    """Print export commands for setting environment variables."""
    env = []
    for name in list_installed_genomes(genomes_dir):
        try:
            g = Genome(name, genomes_dir, build_index=False)
            env_name = re.sub(r"[^\w]+", "_", name).upper()
            env.append(f"export {env_name}={g.filename}")
        except (
            pyfaidx.FastaIndexingError,
            pyfaidx.IndexNotFoundError,
            FileNotFoundError,
        ):
            pass
    return env


def _lazy_provider_selection(name, provider=None):
    """return the first PROVIDER which has genome NAME"""
    providers = []
    for p in online_providers(provider):
        providers.append(p.name)
        if name in p.genomes or (
            p.name == "URL" and try_except_pass(ValueError, check_url, name)
        ):
            return p

    raise GenomeDownloadError(f"{name} not found on {', '.join(providers)}.")


def _provider_selection(name, localname, genomes_dir, provider=None):
    """
    Return a provider object

    First tries to return a specified provider,
    Second tries to return the provider from the README
    Third tries to return the first provider which has the genome (Ensembl>UCSC>NCBI)
    """
    readme = os.path.join(genomes_dir, localname, "README.txt")
    if provider is None and os.path.exists(readme):
        m, _ = read_readme(readme)
        p = m["provider"].lower()
        if p in ["ensembl", "ucsc", "ncbi"]:
            provider = p

    return _lazy_provider_selection(name, provider)


def _filter_genome(
    genome_file: str,
    regex: Optional[str] = None,
    invert_match: Optional[bool] = False,
    keep_alt: Optional[bool] = False,
):
    """
    Combine regex filters & document filtered contigs

    keep_alt : bool , optional
        Set to true to keep these alternative regions.

    regex : str , optional
        Regular expression to select specific chromosome / scaffold names.

    invert_match : bool , optional
        Set to True to select all chromosomes that don't match the regex.
    """
    if keep_alt is True and regex is None:
        return

    regex_func = _get_fasta_regex_func(regex, invert_match, keep_alt)
    excluded_contigs = _apply_fasta_regex_func(genome_file, regex_func)

    # document
    regex_line = "regex: "
    if keep_alt is False:
        regex_line += "'alt' (inverted match)" + (" and " if regex else "")
    if regex:
        regex_line += f"'{regex}'" + (" (inverted match)" if invert_match else "")
    lines = ["", regex_line, ""] + (
        [
            "The following contigs were filtered out of the genome:",
            f"{', '.join(excluded_contigs)}",
        ]
        if excluded_contigs
        else ["No contigs were removed."]
    )

    readme = os.path.join(os.path.dirname(genome_file), "README.txt")
    update_readme(readme, extra_lines=lines)


def _get_fasta_regex_func(
    regex: Optional[str] = None,
    invert_match: Optional[bool] = False,
    keep_alt: Optional[bool] = False,
):
    """
    returns a regex function that accepts a contig header and returns a bool to keep the contig or not.
    """
    # define filter functions
    if keep_alt is False:
        alt_pattern = re.compile("(alt)", re.I)  # case insensitive

        def alt_keep(header):
            return not bool(alt_pattern.search(header))

        keep = alt_keep  # rename in case there is only 1 filter function

    if regex is not None:
        re_pattern = re.compile(regex)

        def re_keep(header):
            return bool(re_pattern.search(header)) is not invert_match

        keep = re_keep  # rename in case there is only 1 filter function

    # combine filter functions?
    if regex is not None and keep_alt is False:

        def keep(header):
            return re_keep(header) and alt_keep(header)

    return keep  # noqa: IDE is confused


def _delete_extensions(directory: str, exts: list):
    """remove (gzipped) files in a directory matching any given extension"""
    for ext in exts:
        [rm_rf(f) for f in glob_ext_files(directory, ext)]


def _is_genome_dir(dirname):
    """
    Check if a directory contains a fasta file of the same name

    Parameters
    ----------
    dirname : str
        Directory name

    Returns
    ------
    bool
    """
    genome_file = os.path.join(dirname, f"{os.path.basename(dirname)}.fa")
    return os.path.exists(genome_file) or os.path.exists(f"{genome_file}.gz")


# Regular expression to check for region (chr:start-end or genome@chr:start-end)
region_p = re.compile(r"^[^@]+@([^\s]+):(\d+)-(\d+)$")


def _check_minsize(fa, minsize):
    """
    Raise ValueError if there is any sequence that is shorter than minsize.
    If minsize is None the size will not be checked.
    """
    if minsize is None:
        return fa

    for name, seq in fa.items():
        if len(seq) < minsize:
            raise ValueError(f"sequence {name} is shorter than {minsize}")

    return fa


def _genomepy_convert(to_convert, genome, minsize=None):
    """
    Convert a variety of inputs using track2fasta().
    """
    if genome is None:
        raise ValueError("input file is not a FASTA file, need a genome!")

    g = Genome(genome)
    tmpfile = NamedTemporaryFile()
    g.track2fasta(to_convert, tmpfile.name)

    fa = as_seqdict(tmpfile.name)
    return _check_minsize(fa, minsize)


def _as_seqdict_genome_regions(regions, minsize=None):
    """
    Accepts list of regions where the genome is encoded in the region,
    using the genome@chrom:start-end format.
    """
    genomic_regions = {}
    for region in regions:
        genome, region = region.split("@")
        genomic_regions.setdefault(genome, []).append(region)

    # test if all genomes are installed
    for genome in genomic_regions:
        Genome(genome)

    tmpfa = NamedTemporaryFile(mode="w", delete=False)
    for genome, g_regions in genomic_regions.items():
        g = Genome(genome)

        fa = g.track2fasta(g_regions)

        for seq in fa:
            seq.name = f"{genome}@{seq.name}"
            print(seq.__repr__(), file=tmpfa)

    tmpfa.flush()

    # Open tempfile and restore original sequence order
    fa = as_seqdict(tmpfa.name)
    fa = {region: fa[region] for region in regions}
    return _check_minsize(fa, minsize)


@singledispatch
def as_seqdict(to_convert, genome=None, minsize=None):
    """
    Convert input to a dictionary with name as key and sequence as value.
    If the input contains genomic coordinates, the genome needs to be
    specified. If minsize is specified all sequences will be checked if they
    are not shorter than minsize. If regions (or a region file) are used as
    the input, the genome can optionally be specified in the region using the
    following format: genome@chrom:start-end.
    Current supported input types include:
    * FASTA, BED and region files.
    * List or numpy.ndarray of regions.
    * pyfaidx.Fasta object.
    * pybedtools.BedTool object.
    Parameters
    ----------
    to_convert : list, str, pyfaidx.Fasta or pybedtools.BedTool
        Input to convert to FASTA-like dictionary
    genome : str, optional
        Genomepy genome name.
    minsize : int or None, optional
        If specified, check if all sequences have at least size minsize.
    Returns
    -------
        dict with sequence names as key and sequences as value.
    """
    raise NotImplementedError(f"Not implement for {type(to_convert)}")


@as_seqdict.register(list)
def _as_seqdict_list(to_convert, genome=None, minsize=None):
    """
    Accepts list of regions as input.
    """
    if region_p.match(to_convert[0]):
        return _as_seqdict_genome_regions(to_convert, minsize)

    return _genomepy_convert(to_convert, genome, minsize)


@as_seqdict.register(TextIOWrapper)
def _as_seqdict_file_object(to_convert, genome=None, minsize=None):
    """
    Accepts file object as input, should be a FASTA file.
    """
    fa = {x: y for x, y in SimpleFastaParser(to_convert)}
    return _check_minsize(fa, minsize)


@as_seqdict.register(str)
def _as_seqdict_filename(to_convert, genome=None, minsize=None):
    """
    Accepts filename as input.
    """
    if not os.path.exists(to_convert):
        raise ValueError("Assuming filename, but it does not exist")

    f = open(to_convert)
    fa = as_seqdict(f)

    if any(fa):
        return _check_minsize(fa, minsize)

    with open(to_convert) as f:
        line = ""
        while True:
            line = f.readline()
            if line == "":
                break
            if not line.startswith("#"):
                break

        if line == "":
            raise IOError(f"empty file {to_convert}")

        if region_p.match(line.strip()):
            regions = [region.strip() for region in [line] + f.readlines()]
            return _as_seqdict_genome_regions(regions, minsize=None)

    # Biopython parser resulted in empty dict
    # Assuming it's a BED or region file
    return _genomepy_convert(to_convert, genome, minsize)


@as_seqdict.register(pyfaidx.Fasta)
def _as_seqdict_pyfaidx(to_convert, genome=None, minsize=None):
    """
    Accepts pyfaidx.Fasta object as input.
    """
    fa = {k: str(v) for k, v in to_convert.items()}
    return _check_minsize(fa, minsize)


try:
    import pybedtools

    @as_seqdict.register(pybedtools.BedTool)
    def _as_seqdict_bedtool(to_convert, genome=None, minsize=None):
        """
        Accepts pybedtools.BedTool as input.
        """
        return _genomepy_convert(
            ["{}:{}-{}".format(*f[:3]) for f in to_convert], genome, minsize
        )


except ImportError:
    pass

try:
    import numpy as np

    @as_seqdict.register(np.ndarray)
    def _as_seqdict_array(to_convert, genome=None, minsize=None):
        """
        Accepts numpy.ndarray with regions as input.
        """
        return as_seqdict(list(to_convert), genome, minsize)


except ImportError:
    pass
