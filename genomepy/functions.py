"""Module-level functions."""
import os
import norns
import re
import sys
from appdirs import user_config_dir, user_cache_dir
from genomepy.__about__ import __version__
from genomepy.exceptions import GenomeDownloadError
from genomepy.genome import Genome
from genomepy.plugin import get_active_plugins, init_plugins
from genomepy.provider import ProviderBase
from genomepy.utils import (
    get_localname,
    get_genomes_dir,
    generate_gap_bed,
    generate_fa_sizes,
    glob_ext_files,
    mkdir_p,
    rm_rf,
    read_readme,
    sanitize_annotation,
    safe,
    check_url,
    try_except_pass,
)
from pyfaidx import FastaIndexingError

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
        except (FastaIndexingError, FileNotFoundError):
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
    keep_alt=False,
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

    keep_alt : bool , optional
        Some genomes contain alternative regions. These regions cause issues with
        sequence alignment, as they are inherently duplications of the consensus regions.
        Set to true to keep these alternative regions.

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
            keep_alt=keep_alt,
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
        if force:
            # overwrite previous versions
            generate_fa_sizes(g.genome_file, g.sizes_file)
            generate_gap_bed(g.genome_file, g.gaps_file)

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
