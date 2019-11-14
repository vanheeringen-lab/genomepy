#!/usr/bin/env python
import click
import genomepy

from collections import deque


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(genomepy.__about__.__version__)
def cli():
    """ Genomes for Python (and others)!

    Version: {}""".format(
        genomepy.__version__
    )
    pass


@click.command("search", short_help="search for genomes")
@click.argument("term")
@click.option("-p", "--provider", help="provider")
def search(term, provider=None):
    """Search for genomes that contain TERM in their name or description."""
    for row in genomepy.search(term, provider):
        print("\t".join([x.decode("utf-8", "ignore") for x in row]))


general_install_options = {
    "genome_dir": {
        "short": "g",
        "long": "genome_dir",
        "help": "genome directory",
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
        "help": "hard/soft/no mask (default: soft)",
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
    "annotation": {
        "short": "a",
        "long": "annotation",
        "help": "download annotation",
        "flag_value": True,
    },
    "force": {
        "short": "f",
        "long": "force",
        "help": "overwrite existing files",
        "flag_value": True,
    },
}


def get_install_options():
    """combine general and provider specific options

    add provider in front of the provider specific options to prevent overlap"""
    install_options = general_install_options

    for name in genomepy.provider.ProviderBase.list_providers():
        p = genomepy.provider.ProviderBase.create(name)
        p_dict = p.list_install_options()
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
    genome_dir,
    localname,
    mask,
    regex,
    invert_match,
    bgzip,
    annotation,
    force,
    **kwargs
):
    """Install genome NAME from provider PROVIDER in directory GENOME_DIR."""
    genomepy.install_genome(
        name,
        provider,
        genome_dir=genome_dir,
        localname=localname,
        mask=mask,
        regex=regex,
        invert_match=invert_match,
        bgzip=bgzip,
        annotation=annotation,
        force=force,
        **kwargs
    )


@click.command("genomes", short_help="list available genomes")
@click.option("-p", "--provider", help="provider")
def genomes(provider=None):
    """List all available genomes."""
    for row in genomepy.list_available_genomes(provider):
        print("\t".join(row))


@click.command("providers", short_help="list available providers")
def providers():
    """List all available providers."""
    for p in genomepy.list_available_providers():
        print(p)


@click.command("plugin", short_help="manage plugins")
@click.argument("command")
@click.argument("name", nargs=-1)
def plugin(command, name):
    """Enable or disable plugins

    Use 'genomepy plugin list' to show all available plugins

    Use 'genomepy plugin enable/disable [NAME]' to (dis)able plugins"""
    genomepy.functions.manage_plugins(command, name)


@click.command("config", short_help="manage configuration")
@click.argument("command")
def config(command):
    """Manage configuration"""
    genomepy.functions.manage_config(command)


cli.add_command(search)
cli.add_command(install)
cli.add_command(genomes)
cli.add_command(providers)
cli.add_command(plugin)
cli.add_command(config)

if __name__ == "__main__":
    cli()
