#!/usr/bin/env python
import os
import sys
from collections import deque

import click
from colorama import Fore, Style, init
from loguru import logger

import genomepy

init(autoreset=True)

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(genomepy.__version__)
def cli():
    """ Genomes for Python (and others)!

    Version: {}""".format(
        genomepy.__version__
    )
    pass  # noqa


@click.command("clean", short_help="remove provider data")
def clean():
    """
    Remove cached data on providers (e.g. available genomes).
    """
    genomepy.clean()


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
    """
    List all available genomes.

    Returns the metadata of each found genome, including the availability of a gene annotation.
    For UCSC, up to 4 gene annotation styles may be available: UCSC, Ensembl, NCBI_refseq and UCSC_refseq.
    """
    terminal_header()
    for row in genomepy.list_available_genomes(provider):
        terminal_formatting(row)


# extended options for genomepy install
INSTALL_OPTIONS = {
    "provider": {
        "short": "p",
        "long": "provider",
        "help": "download from this provider",
        "default": None,
    },
    "genomes_dir": {
        "short": "g",
        "long": "genomes_dir",
        "help": "create output directory here",
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
    "keep_alt": {
        "short": "k",
        "long": "keep-alt",
        "help": "keep alternative regions",
        "flag_value": True,
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
    "skip_matching": {
        "short": "sm",
        "long": "skip_matching",
        "help": "skip matching contigs between the gene annotation and the genome (sets -a)",
        "flag_value": True,
    },
    "skip_filter": {
        "short": "sf",
        "long": "skip_filter",
        "help": "skip filtering out contigs in the gene annotation missing from the genome (sets -a)",
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
    install_options = INSTALL_OPTIONS

    # extend install options with provider specific options
    if "install" in click.get_os_args():
        for provider in genomepy.list_providers():
            p_dict = eval(
                "genomepy.provider."
                + provider.capitalize()
                + "Provider.provider_specific_install_options"
            )
            for option in p_dict.keys():
                p_dict[option]["long"] = provider + "-" + p_dict[option]["long"]
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
@click.argument("name")
@cli.command()
def install(
    name,
    provider,
    genomes_dir,
    localname,
    mask,
    keep_alt,
    regex,
    invert_match,
    bgzip,
    annotation,
    only_annotation,
    skip_matching,
    skip_filter,
    threads,
    force,
    **kwargs,
):
    """install a genome & run active plugins"""
    genomepy.install_genome(
        name,
        provider=provider,
        genomes_dir=genomes_dir,
        localname=localname,
        mask=mask,
        keep_alt=keep_alt,
        regex=regex,
        invert_match=invert_match,
        bgzip=bgzip,
        annotation=annotation,
        only_annotation=only_annotation,
        skip_matching=skip_matching,
        skip_filter=skip_filter,
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
    for p in genomepy.list_providers():
        print(p)


# names and sizes for the search output columns (when connected to a terminal)
SEARCH_FORMAT = {
    "name": 20,
    "provider": 8,  # fixed width
    "accession": 16,  # fixed width
    "tax_id": 7,  # fixed width
    "annotation": 10,  # fixed width
    "species": 40,
    "other_info": 40,
}
SEARCH_STRING = "    ".join([f"{{: <{size}}}" for size in SEARCH_FORMAT.values()])
if sys.stdout.isatty():

    def bool_to_unicode(boolean: bool) -> str:
        """converts True to a checkmark and False to a cross-mark"""
        return "\u2713" if boolean else "\u2717"

    def color_unicode(string):
        sting = string.replace("\u2713", Fore.GREEN + "\u2713" + Fore.RESET)
        sting = sting.replace("\u2717", Fore.RED + "\u2717" + Fore.RESET)
        return sting

    def terminal_formatting(row: list):
        """
        In case we print to a terminal, the output is aligned.
        Otherwise (file, pipe) we use tab-separated columns.
        """
        if isinstance(row[4], list):
            row[4] = " ".join([bool_to_unicode(b) for b in row[4]])
        else:
            row[4] = bool_to_unicode(row[4])
        for n, ele in enumerate(row):
            if ele is None:
                row[n] = "na"
        row = SEARCH_STRING.format(*row)
        print(color_unicode(row))

    def terminal_header():
        print(Style.BRIGHT + SEARCH_STRING.format(*SEARCH_FORMAT))


else:

    def terminal_formatting(row: list):
        """
        In case we print to a terminal, the output is aligned.
        Otherwise (file, pipe) we use tab-separated columns.
        """
        if isinstance(row[4], list):
            row[4] = str(row[4])
        print("\t".join([str(element) for element in row]))

    def terminal_header():
        print("\t".join(SEARCH_FORMAT))


@click.command("search", short_help="search for genomes")
@click.argument("term")
@click.option("-p", "--provider", help="provider")
def search(term, provider=None):
    """
    Search for genomes that contain TERM in their name, description,
    accession (must start with GCA_ or GCF_) or taxonomy (exact matches only).

    Function is case-insensitive. Spaces in TERM can be replaced with underscores
    (_) or TERM can be "quoted", e.g., "homo sapiens".

    Returns the metadata of each found genome, including the availability of a gene annotation.
    For UCSC, up to 4 gene annotation styles may be available: UCSC, Ensembl, NCBI_refseq and UCSC_refseq.
    """
    no_genomes = True
    for row in genomepy.search(term, provider):
        if no_genomes:
            no_genomes = False
            terminal_header()
        terminal_formatting(row)

    if sys.stdout.isatty():
        if no_genomes:
            logger.error("No genomes found!")
        else:
            print(Fore.GREEN + " ^")
            print(Fore.GREEN + " Use name for " + Fore.CYAN + "genomepy install")


cli.add_command(clean)
cli.add_command(config)
cli.add_command(genomes)
cli.add_command(install)
cli.add_command(plugin)
cli.add_command(providers)
cli.add_command(search)

if __name__ == "__main__":
    cli()
