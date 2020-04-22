#!/usr/bin/env python
import click
import genomepy
import sys
import os

from collections import deque
from colorama import init, Fore, Style

init(autoreset=True)

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(genomepy.__about__.__version__)
def cli():
    """ Genomes for Python (and others)!

    Version: {}""".format(
        genomepy.__version__
    )
    pass


@click.command("config", short_help="manage configuration")
@click.argument("command")
def config(command):
    """
    Manage configuration

    genomepy config file        return config filepath

    genomepy config show        return config content

    genomepy config generate    create new config file
    """
    genomepy.manage_config(command)


@click.command("genomes", short_help="list available genomes")
@click.option("-p", "--provider", help="provider")
def genomes(provider=None):
    """List all available genomes."""
    for row in genomepy.list_available_genomes(provider):
        print("\t".join(row))


# extended options for genomepy install
general_install_options = {
    "genomes_dir": {
        "short": "g",
        "long": "genomes_dir",
        "help": "genomes directory",
        "default": None,
    },
    "localname": {
        "short": "l",
        "long": "localname",
        "help": "custom name",
        "default": None,
    },
    "mask": {
        "short": "m",
        "long": "mask",
        "help": "DNA masking: hard/soft/none (default: soft)",
        "default": "soft",
    },
    "regex": {
        "short": "r",
        "long": "regex",
        "help": "regex to filter sequences",
        "default": None,
    },
    "invert_match": {
        "short": "n",
        "long": "no-match",
        "help": "select sequences that *don't* match regex",
        "flag_value": True,
    },
    "bgzip": {
        "short": "b",
        "long": "bgzip",
        "help": "bgzip genome",
        "flag_value": True,
    },
    "threads": {
        "short": "t",
        "long": "threads",
        "help": "build index using multithreading",
        "default": min(os.cpu_count(), 8),
    },
    "force": {
        "short": "f",
        "long": "force",
        "help": "overwrite existing files",
        "flag_value": True,
    },
    "text_line1": {
        "long": "Annotation options:",
        "help": "",
        "flag_value": True,
        "text_line": True,
    },
    "annotation": {
        "short": "a",
        "long": "annotation",
        "help": "download annotation",
        "flag_value": True,
    },
    "only_annotation": {
        "short": "o",
        "long": "only_annotation",
        "help": "only download annotation (sets -a)",
        "flag_value": True,
    },
    "skip_sanitizing": {
        "short": "s",
        "long": "skip_sanitizing",
        "help": "skip (check for) matching of contig names between annotation and fasta (sets -a)",
        "flag_value": True,
    },
    "text_line2": {
        "long": "Provider specific options:",
        "help": "",
        "flag_value": True,
        "text_line": True,
    },
}


def get_install_options():
    """combine general and provider specific options

    add provider in front of the provider specific options to prevent overlap"""
    install_options = general_install_options

    for name in genomepy.ProviderBase.list_providers():
        p_dict = genomepy.ProviderBase.create(name).list_install_options()
        for option in p_dict.keys():
            p_dict[option]["long"] = name + "-" + p_dict[option]["long"]
        install_options.update(p_dict)

    return install_options


def custom_options(options):
    """dynamically add options to a click.command (based on a dict with options)"""

    def decorator(f):
        for opt_name, opt_params in options.items():
            param_decls = deque(["--" + opt_params["long"], opt_name])
            if "short" in opt_params.keys():
                param_decls.appendleft("-" + opt_params["short"])

            attrs = dict(help=opt_params["help"])
            if "type" in opt_params.keys():
                attrs["type"] = opt_params["type"]
            if "default" in opt_params.keys():
                attrs["default"] = opt_params["default"]
            if "flag_value" in opt_params.keys():
                attrs["flag_value"] = opt_params["flag_value"]

            # can be used to add paragraphs to the --help menu
            if "text_line" in opt_params.keys():
                param_decls = deque(["\n" + opt_params["long"], opt_name])

            click.option(*param_decls, **attrs)(f)
        return f

    return decorator


@custom_options(get_install_options())
@click.argument("provider")
@click.argument("name")
@cli.command()
def install(
    name,
    provider,
    genomes_dir,
    localname,
    mask,
    regex,
    invert_match,
    bgzip,
    annotation,
    only_annotation,
    skip_sanitizing,
    threads,
    force,
    **kwargs,
):
    """Install genome NAME from provider PROVIDER in directory GENOME_DIR."""
    genomepy.install_genome(
        name,
        provider,
        genomes_dir=genomes_dir,
        localname=localname,
        mask=mask,
        regex=regex,
        invert_match=invert_match,
        bgzip=bgzip,
        annotation=annotation,
        only_annotation=only_annotation,
        skip_sanitizing=skip_sanitizing,
        threads=threads,
        force=force,
        **kwargs,
    )


@click.command("plugin", short_help="manage plugins")
@click.argument("command")
@click.argument("name", nargs=-1)
def plugin(command, name):
    """
    Enable or disable plugins

    genomepy plugin list                 show plugins and status

    genomepy plugin enable  [NAME(S)]    enable plugins

    genomepy plugin disable [NAME(S)]    disable plugins
    """
    genomepy.manage_plugins(command, name)


@click.command("providers", short_help="list available providers")
def providers():
    """List all available providers."""
    for p in genomepy.list_available_providers():
        print(p)


@click.command("search", short_help="search for genomes")
@click.argument("term")
@click.option("-p", "--provider", help="provider")
def search(term, provider=None):
    """
    Search for genomes that contain TERM in their name or description.

    Function is case-insensitive. Spaces in TERM can be replaced with underscores
    (_) or TERM can be "quoted", e.g., "homo sapiens".
    """
    data = [["name", "provider", "accession", "species", "tax_id", "other_info"]]
    for row in genomepy.search(term, provider):
        data.append([x.decode("utf-8", "ignore") for x in row])
    if len(data) == 1:
        print("No genomes found!", file=sys.stderr)
        return

    # In case we print to a terminal, the output is aligned.
    # Otherwise (file, pipe) we use tab-separated columns.
    if sys.stdout.isatty():
        sizes = [max(len(row[i]) + 4 for row in data) for i in range(len(data[0]))]
        fstring = "".join([f"{{: <{size}}}" for size in sizes])
    else:
        fstring = "\t".join(["{}" for _ in range(len(data[0]))])

    for i, row in enumerate(data):
        if i == 0:
            print(Style.BRIGHT + fstring.format(*row))
        else:
            print(fstring.format(*row))
    if sys.stdout.isatty():
        print(Fore.GREEN + " ^")
        print(Fore.GREEN + " Use name for " + Fore.CYAN + "genomepy install")


cli.add_command(config)
cli.add_command(genomes)
cli.add_command(install)
cli.add_command(plugin)
cli.add_command(providers)
cli.add_command(search)

if __name__ == "__main__":
    cli()
