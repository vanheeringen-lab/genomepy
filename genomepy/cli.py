#!/usr/bin/env python
import click
import os
import sys
import genomepy

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(genomepy.__about__.__version__)
def cli():
    """ Genomes for Python (and others)!
    
    Version: {}""".format(genomepy.__version__)
    pass

@click.command('search', short_help="search for genomes")
@click.argument("term")
@click.option("-p", "--provider", help="provider")
def search(term, provider=None):
    """Search for genomes that contain TERM in their name or description."""
    for row in genomepy.search(term, provider):
        print("\t".join([x.decode('utf-8', 'ignore') for x in row]))

@click.command('install', short_help="install genome")
@click.argument("name")
@click.argument("provider")
@click.option("-g", "--genome_dir", help="genome directory", default=None)
@click.option("-l", "--localname", help="custom name", default=None)
@click.option("-m", "--mask", help="mask (hard or soft)", default="soft")
@click.option("-r", "--regex", help="regex to filter sequences", default=None)
@click.option("--match/--no-match", help="set no-match to select sequences that *don't*  match regex", default=True)
@click.option("--annotation/--no-annotation", help="download annotation", default=False)
def install(name, provider, genome_dir, localname, mask, regex, match, annotation):
    """Install genome NAME from provider PROVIDER in directory GENOME_DIR."""
    genomepy.install_genome(
            name, provider, genome_dir=genome_dir, localname=localname, mask=mask, 
            regex=regex, invert_match=not(match), annotation=annotation)

@click.command('genomes', short_help="list available genomes")
@click.option("-p", "--provider", help="provider")
def genomes(provider=None):
    """List all available genomes."""
    for row in genomepy.list_available_genomes(provider):
        print("\t".join(row))

@click.command('providers', short_help="list available providers")
def providers():
    """List all available providers."""
    for p in genomepy.list_available_providers():
        print(p)

@click.command('plugin', short_help="manage plugins")
@click.argument("command")
@click.argument("name", nargs=-1)
def plugin(command, name):
    """Enable or disable plugins"""
    genomepy.functions.manage_plugins(command, name)

@click.command('config', short_help="manage configuration")
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
