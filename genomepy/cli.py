#!/usr/bin/env python
import click
import genomepy


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


install_options = {
    "genome_dir": {
        "short": "g",
        "long": "genome_dir",
        "help": "genome directory",
        # "type": str,
        "default": None,
    },
    "localname": {
        "short": "l",
        "long": "localname",
        "help": "custom name",
        # "type": str,
        "default": None,
    },
    "mask": {
        "short": "m",
        "long": "mask",
        "help": "mask (hard/soft/none) default=soft",
        # "type": str,
        "default": "soft",
    },
    "regex": {
        "short": "r",
        "long": "regex",
        "help": "regex to filter sequences",
        # "type": str,
        "default": None,
    },
    "invert_match": {
        "short": "n",
        "long": "no-match",
        "help": "select sequences that *don't* match regex",
        # "type": bool,
        "flag_value": True,
    },
    "annotation": {
        "short": "a",
        "long": "annotation",
        "help": "download annotation",
        # "type": bool,
        "flag_value": True,
    },
    "force": {
        "short": "f",
        "long": "force",
        "help": "overwrite existing files",
        # "type": bool,
        "flag_value": True,
    },
    "toplevel": {
        "short": "t",
        "long": "toplevel",
        "help": "always download toplevel-genome (Ensembl)",
        # "type": bool,
        "flag_value": True,
    },
}


def custom_options(options):
    def decorator(f):
        for opt_name, opt_params in options.items():
            param_decls = (
                "-" + opt_params["short"],
                "--" + opt_params["long"],
                opt_name,  # opt_params["name"],
            )
            attrs = dict(help=opt_params["help"])
            # type=opt_params["type"])
            if "type" in opt_params.keys():
                attrs["type"] = opt_params["type"]
            if "default" in opt_params.keys():
                attrs["default"] = opt_params["default"]
            if "flag_value" in opt_params.keys():
                attrs["flag_value"] = opt_params["flag_value"]

            click.option(*param_decls, **attrs)(f)
        return f

    return decorator


# provider_options = {
#     'Ensembl': {
#         'toplevel': ["BOOL", "always download toplevel-genome", "flag_value=True"],
#         'version': ["INT", "select release version", "None"]},
#     'NCBI': {},
#     'UCSC': {}}
# for provider, options in provider_options.items():
#     for option, info in options.items():
#         print("  --{}-{} {}\t{}".format(provider, option, info[0], info[1]))
#         # @click.option("--{}-{}".format(provider, option), help=info[1], eval(info[2]))


# def get_PSOs():
#     PSOs = ['Ensembl-version=VERSION', 'Ensembl-toplevel']
#     parsed_list = "\n".join(PSOs)
#     return parsed_list


@click.group()
def cli():
    pass


@custom_options(install_options)
@click.argument("name")
@click.argument("provider")
@cli.command()
# @click.command("install", short_help="install genome")
# @click.argument("name")
# @click.argument("provider")
# # @click.argument("provider_specific_options", nargs=-1)
# @click.option("-g", "--genome_dir", help="genome directory", default=None)
# @click.option("-l", "--localname", help="custom name", default=None)
# @click.option(
#     "-m", "--mask", help="mask (hard, soft or none. Default = soft)", default="soft"
# )
# @click.option("-r", "--regex", help="regex to filter sequences", default=None)
# @click.option(
#     "-n",
#     "--no_match",
#     "invert_match",
#     help="select sequences that *don't* match regex",
#     flag_value=True,
# )
# @click.option(
#     "-t",
#     "--toplevel",
#     help="always download toplevel-genome (Ensembl)",
#     flag_value=True,
# )
# @click.option("-a", "--annotation", help="download annotation", flag_value=True)
# @click.option("-f", "--force", help="overwrite existing files", flag_value=True)
# # @click.option("--PSO", help="Provider specific options:\n{}".format(get_PSOs()), flag_value=True)
def install(
    name,
    provider,
    genome_dir,
    localname,
    mask,
    toplevel,
    regex,
    force,
    invert_match,
    annotation,
):
    """Install genome NAME from provider PROVIDER in directory GENOME_DIR."""

    # """Install genome NAME from provider PROVIDER in directory GENOME_DIR.
    #
    # PSOs:\n{}""".format(get_PSOs())

    # """Install genome NAME from provider PROVIDER in directory GENOME_DIR.
    #
    # PSOs:"""+get_PSOs()
    genomepy.install_genome(
        name,
        provider,
        genome_dir=genome_dir,
        localname=localname,
        mask=mask,
        toplevel=toplevel,
        regex=regex,
        force=force,
        invert_match=invert_match,
        annotation=annotation,
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

    Use 'genomepy enable/disable [NAME]' to (dis)able plugins"""
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
