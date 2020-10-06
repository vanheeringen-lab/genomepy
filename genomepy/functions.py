"""Module-level functions."""
import os
import re
import sys
from functools import singledispatch
from io import TextIOWrapper

from appdirs import user_config_dir, user_cache_dir
import norns
import pyfaidx
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tempfile import NamedTemporaryFile

from genomepy.__about__ import __version__
from genomepy.exceptions import GenomeDownloadError
from genomepy.genome import Genome
from genomepy.plugin import get_active_plugins, init_plugins
from genomepy.provider import ProviderBase
from genomepy.utils import (
    get_localname,
    get_genomes_dir,
    glob_ext_files,
    mkdir_p,
    rm_rf,
    read_readme,
    sanitize_annotation,
    safe,
    check_url,
    try_except_pass,
)

config = norns.config("genomepy", default="cfg/default.yaml")


def clean():
    """Remove cached data on providers"""
    my_cache_dir = os.path.join(user_cache_dir("genomepy"), __version__)
    rm_rf(my_cache_dir)
    mkdir_p(my_cache_dir)
    print("All clean!")


def manage_config(cmd):
    """Manage genomepy config file."""
    if cmd == "file":
        print(config.config_file)
    elif cmd == "show":
        with open(config.config_file) as f:
            print(f.read())
    elif cmd == "generate":
        config_dir = user_config_dir("genomepy")
        if not os.path.exists(config_dir):
            mkdir_p(config_dir)

        new_config = os.path.join(config_dir, "genomepy.yaml")
        # existing config must be removed before norns picks up the default again
        if os.path.exists(new_config):
            os.unlink(new_config)
        default_config = norns.config(
            "genomepy", default="cfg/default.yaml"
        ).config_file
        with open(new_config, "w") as fout, open(default_config) as fin:
            fout.write(fin.read())
        config.config_file = new_config
        print(f"Created config file {new_config}")
    else:
        raise ValueError(f"Invalid config command: {cmd}")


def _online_providers():
    """Return a list of online providers as objects"""
    providers = []
    for p in ProviderBase.list_providers():
        try:
            providers.append(ProviderBase.create(p))
        except ConnectionError as e:
            sys.stderr.write(str(e))
    return providers


def _providers(provider=None):
    """
    Return a list of provider objects:
    either the specified provider, or all online providers
    """
    return [ProviderBase.create(provider)] if provider else _online_providers()


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
    providers = _providers(provider)
    for p in providers:
        for row in p.list_available_genomes():
            yield [p.name] + list(row)


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
    name = os.path.basename(os.path.abspath(os.path.expanduser(dirname)))
    return any([f for f in glob_ext_files(dirname) if os.path.join(name, name) in f])


def list_installed_genomes(genomes_dir=None):
    """
    List all available genomes.

    Parameters
    ----------
    genomes_dir : str
        Directory with installed genomes.

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


def generate_exports(genomes_dir=None):
    """Print export commands for setting environment variables."""
    env = []
    for name in list_installed_genomes(genomes_dir):
        try:
            g = Genome(name)
            env_name = re.sub(r"[^\w]+", "_", name).upper()
            env.append(f"export {env_name}={g.filename}")
        except (pyfaidx.FastaIndexingError, FileNotFoundError):
            pass
    return env


def generate_env(fname="exports.txt", genomes_dir=None):
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
        for env in generate_exports(genomes_dir):
            fout.write(f"{env}\n")


def _lazy_provider_selection(name, provider=None):
    """return the first PROVIDER which has genome NAME"""
    providers = _providers(provider)
    for p in providers:
        if name in p.genomes or (
            p.name == "URL" and try_except_pass(ValueError, check_url, name)
        ):
            return p

    raise GenomeDownloadError(
        f"{name} not found on {', '.join([p.name for p in providers])}."
    )


def _provider_selection(name, localname, genomes_dir, provider=None):
    """
    Return a provider object

    First tries to return a specified provider,
    Second tries to return the provider from the README
    Third tries to return the first provider which has the genome (Ensembl>UCSC>NCBI)
    """
    if provider is None:
        readme = os.path.join(genomes_dir, localname, "README.txt")
        m, _ = read_readme(readme)
        p = m["provider"].lower()
        if p in ["ensembl", "ucsc", "ncbi"]:
            provider = p

    return _lazy_provider_selection(name, provider)


def install_genome(
    name,
    provider=None,
    genomes_dir=None,
    localname=None,
    mask="soft",
    regex=None,
    invert_match=False,
    bgzip=None,
    annotation=False,
    only_annotation=False,
    skip_sanitizing=False,
    threads=1,
    force=False,
    **kwargs,
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

    regex : str , optional
        Regular expression to select specific chromosome / scaffold names.

    invert_match : bool , optional
        Set to True to select all chromosomes that don't match the regex.

    bgzip : bool , optional
        If set to True the genome FASTA file will be compressed using bgzip.
        If not specified, the setting from the configuration file will be used.

    threads : int , optional
        Build genome index using multithreading (if supported). Default: lowest of 8/all threads

    force : bool , optional
        Set to True to overwrite existing files.

    annotation : bool , optional
        If set to True, download gene annotation in BED and GTF format.

    only_annotation : bool , optional
        If set to True, only download the annotation files.

    skip_sanitizing : bool , optional
        If set to True, downloaded annotation files whose sequence names do not match
        with the (first header fields of) the genome.fa will not be corrected.

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

    # Check if genome already exists, or if downloading is forced
    genome_found = _is_genome_dir(out_dir)
    if (not genome_found or force) and not only_annotation:
        # Download genome from provider
        p = _provider_selection(name, localname, genomes_dir, provider)
        p.download_genome(
            name,
            genomes_dir,
            mask=mask,
            regex=regex,
            invert_match=invert_match,
            localname=localname,
            bgzip=bgzip,
            **kwargs,
        )
        genome_found = True

        # Export installed genome(s)
        generate_env(genomes_dir=genomes_dir)

    # Generates a Fasta object, index, gaps and sizes file
    g = None
    if genome_found:
        g = Genome(localname, genomes_dir=genomes_dir)

    # Check if any annotation flags are given, if annotation already exists, or if downloading is forced
    if any(
        [
            annotation,
            only_annotation,
            skip_sanitizing,
            kwargs.get("to_annotation"),
            kwargs.get("ucsc_annotation_type"),
        ]
    ):
        annotation = True
    annotation_found = bool(glob_ext_files(out_dir, "gtf"))
    if (not annotation_found or force) and annotation:
        # Download annotation from provider
        p = _provider_selection(name, localname, genomes_dir, provider)
        p.download_annotation(name, genomes_dir, localname=localname, **kwargs)

        # Sanitize annotation if needed (requires genome)
        annotation_found = bool(glob_ext_files(out_dir, "gtf"))
        if genome_found and annotation_found and not skip_sanitizing:
            sanitize_annotation(g)

    if genome_found:
        # Run all active plugins (requires genome)
        for plugin in get_active_plugins():
            plugin.after_genome_download(g, threads, force)


def manage_plugins(command, plugin_names=None):
    """List, enable or disable plugins."""
    plugins = init_plugins()
    for name in plugin_names if plugin_names else []:
        if name not in plugins:
            raise ValueError(f"Unknown plugin: {name}")

    active_plugins = config.get("plugin", [])
    if command == "list":
        print("{:20}{}".format("plugin", "enabled"))
        for plugin in sorted(plugins):
            print(
                "{:20}{}".format(
                    plugin, {False: "", True: "*"}[plugin in active_plugins]
                )
            )
        return

    elif command in ["enable", "activate"]:
        for name in plugin_names:
            if name not in active_plugins:
                active_plugins.append(name)

    elif command in ["disable", "deactivate"]:
        for name in plugin_names:
            if name in active_plugins:
                active_plugins.remove(name)

    else:
        raise ValueError(f"Invalid plugin command: {command}")

    config["plugin"] = active_plugins
    config.save()
    print("Enabled plugins: {}".format(", ".join(sorted(active_plugins))))


def list_available_providers():
    """
    List all available providers.

    Returns
    -------
    list with provider names
    """
    return ProviderBase.list_providers()


def search(term, provider=None):
    """
    Search for a genome.

    If provider is specified, search only that specific provider, else
    search all providers. Both the name and description are used for the
    search. Search term is case-insensitive.

    Parameters
    ----------
    term : str or int
        Search term, case-insensitive.

    provider : str , optional
        Provider name

    Yields
    ------
    tuple
        genome information (name/identifier and description)
    """
    term = safe(str(term))
    providers = _providers(provider)
    for p in providers:
        for row in p.search(term):
            yield [x.encode("utf-8") for x in list(row[:1]) + [p.name] + list(row[1:])]


# Regular expression to check for region (chr:start-end or genome@chr:start-end)
region_p = re.compile(r"^[^@]+@([^\s]+):(\d+)-(\d+)$")


def _check_minsize(fa, minsize):
    """
    Raise ValueError if there is any sequence that is shorter than minsize.
    If minsize is None the size will not be checked.
    """
    if minsize is not None:
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
        if genome not in genomic_regions:
            Genome(genome)
            genomic_regions[genome] = []
        genomic_regions[genome].append(region)

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
