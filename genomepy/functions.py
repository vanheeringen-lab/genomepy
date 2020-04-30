"""Module-level functions."""
import os
import norns
import re

from appdirs import user_config_dir
from glob import glob
from pyfaidx import FastaIndexingError
from genomepy.genome import Genome
from genomepy.provider import ProviderBase
from genomepy.plugin import get_active_plugins, init_plugins
from genomepy.utils import (
    get_localname,
    get_genomes_dir,
    glob_ext_files,
    mkdir_p,
    sanitize_annotation,
)

config = norns.config("genomepy", default="cfg/default.yaml")


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
    if provider:
        providers = [ProviderBase.create(provider)]
    else:
        # if provider is not specified search all providers
        providers = [ProviderBase.create(p) for p in ProviderBase.list_providers()]

    for p in providers:
        for row in p.list_available_genomes():
            yield [p.name] + list(row)


def _is_genome_dir(dirname):
    """
    Check if a directory contains a fasta file

    Parameters
    ----------
    dirname : str
        Directory name

    Returns
    ------
    bool
    """
    return len(glob(f"{dirname}/*.fa")) > 0


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

    return (
        [
            f
            for f in os.listdir(genomes_dir)
            if _is_genome_dir(os.path.join(genomes_dir, f))
        ]
        if os.path.exists(genomes_dir)
        else []
    )


def generate_exports():
    """Print export commands for setting environment variables."""
    env = []
    for name in list_installed_genomes():
        try:
            g = Genome(name)
            env_name = re.sub(r"[^\w]+", "_", name).upper()
            env.append(f"export {env_name}={g.filename}")
        except FastaIndexingError:
            pass
    return env


def generate_env(fname=None):
    """Generate file with exports.

    By default this is .config/genomepy/exports.txt.

    An alternative file name or file path is accepted too.

    Parameters
    ----------
    fname: str, optional
        Absolute path or name of the output file.
    """
    path_name = os.path.expanduser(str(fname))
    if fname and os.path.exists(os.path.dirname(path_name)):
        absname = os.path.abspath(path_name)
    else:
        config_dir = user_config_dir("genomepy")
        if not os.path.exists(config_dir):
            manage_config("generate")

        name = "exports.txt" if fname is None else fname
        absname = os.path.join(config_dir, name)

    with open(absname, "w") as fout:
        for env in generate_exports():
            fout.write(f"{env}\n")


def install_genome(
    name,
    provider,
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

    provider : str
        Provider name

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
    genomes_dir = get_genomes_dir(genomes_dir, check_exist=False)
    localname = get_localname(name, localname)
    out_dir = os.path.join(genomes_dir, localname)

    # Check if genome already exists, or if downloading is forced
    genome_found = (
        len([f for f in glob_ext_files(out_dir) if f"{localname}.fa" in f]) >= 1
    )
    if (not genome_found or force) and not only_annotation:
        # Download genome from provider
        p = ProviderBase.create(provider)
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
        generate_env()

    # Generates a Fasta object, index, gaps and sizes file
    g = None
    if genome_found:
        g = Genome(localname, genomes_dir=genomes_dir)

    # Check if any annotation flags are given, if annotation already exists, or if downloading is forced
    if any([annotation, only_annotation, kwargs.get("to_annotation", False)]):
        annotation = True
    annotation_found = len(glob_ext_files(out_dir, "gtf")) >= 1
    if (not annotation_found or force) and annotation:
        # Download annotation from provider
        p = ProviderBase.create(provider)
        p.download_annotation(name, genomes_dir, localname=localname, **kwargs)

        # Sanitize annotation if needed (requires genome)
        annotation_found = len(glob_ext_files(out_dir, "gtf")) >= 1
        if genome_found and annotation_found and not skip_sanitizing:
            sanitize_annotation(g)

    if genome_found:
        # Run all active plugins (requires genome)
        for plugin in get_active_plugins():
            plugin.after_genome_download(g, threads, force)


def manage_plugins(command, plugin_names=None):
    """List, enable or disable plugins.
    """
    if command not in ["list", "enable", "disable"]:
        raise ValueError(f"Invalid plugin command: {command}")

    plugins = init_plugins()
    active_plugins = config.get("plugin", [])

    if command == "list":
        print("{:20}{}".format("plugin", "enabled"))
        for plugin in sorted(plugins):
            print(
                "{:20}{}".format(
                    plugin, {False: "", True: "*"}[plugin in active_plugins]
                )
            )
    else:
        if plugin_names:
            for name in plugin_names:
                if name not in plugins:
                    raise ValueError(f"Unknown plugin: {name}")
        else:
            plugin_names = []

        if command == "enable":
            for name in plugin_names:
                if name not in active_plugins:
                    active_plugins.append(name)
        elif command == "disable":
            for name in plugin_names:
                if name in active_plugins:
                    active_plugins.remove(name)

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
    term : str
        Search term, case-insensitive.

    provider : str , optional
        Provider name

    Yields
    ------
    tuple
        genome information (name/identfier and description)
    """
    if provider:
        providers = [ProviderBase.create(provider)]
    else:
        # if provider is not specified search all providers
        providers = [ProviderBase.create(p) for p in ProviderBase.list_providers()]
    for p in providers:
        for row in p.search(term):
            yield [
                x.encode("latin-1") for x in list(row[:1]) + [p.name] + list(row[1:])
            ]
