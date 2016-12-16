import os
import glob
from pyfaidx import Fasta
from genomepy.provider import ProviderBase
__version__ = '0.0.2'
__author__ = "Simon van Heeringen"

def list_available_genomes(provider=None):
    """
    List all available genomes.

    Parameters
    ----------
    provider : str, optional
        List genomes from specific provider. Genomes from all
        provider will be returned if not specified.

    Returns
    -------
    list with genome names
    """
    if provider:
        providers = [ProviderBase.create(provider)]
    else:
        # if provider is not specified search all providers
        providers = [ProviderBase.create(p) for 
                        p in ProviderBase.list_providers()]

    for p in providers:
        for row in p.list_available_genomes():
            yield [p._name] + list(row)

def list_available_providers():
    """
    List all available providers.

    Returns
    -------
    list with provider names
    """
    return ProviderBase.list_providers()

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
    return len(glob.glob("{}/*.fa".format(dirname))) > 0

def list_installed_genomes(genome_dir):
    """
    List all available genomes.

    Parameters
    ----------
    genome_dir : str
        Directory with installed genomes.

    Returns
    -------
    list with genome names
    """

    return [f for f in os.listdir(genome_dir) if 
            _is_genome_dir(genome_dir + "/" + f)]

def search(term, provider=None):
    """
    Search for a genome. If provider is specified, search that 
    specific provider, else search all providers."

    Parameters
    ----------
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
        providers = [ProviderBase.create(p) for 
                        p in ProviderBase.list_providers()]

    for p in providers:
        for row in p.search(term):
            yield [p._name] + list(row)

def install_genome(name, provider, genome_dir):
    """
    Install a genome.

    Parameters
    ----------
    name : str
        Genome name

    provider : str
        Provider name

    genome_dir : str
        Where to store the fasta files
    """
    
    p = ProviderBase.create(provider)
    p.download_genome(name, genome_dir)

def genome(name, genome_dir):
    """
    Get pyfaidx Fasta object of genome

    Parameters
    ----------
    name : str
        Genome name
    
    genome_dir : str
        Genome installation directory

    Returns
    -------
    pyfaidx.Fasta object
    """

    fa = glob.glob("{}/{}/*.fa".format(genome_dir, name))[0]
    return Fasta(fa)
