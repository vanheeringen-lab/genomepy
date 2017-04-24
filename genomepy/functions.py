"""Module-level functions."""
import os
import glob
from pyfaidx import Fasta
from genomepy.provider import ProviderBase
from genomepy.utils import generate_sizes
import norns

config = norns.config("genomepy", default="cfg/default.yaml")

# Python 2
try:
    FileNotFoundError
except NameError:
    # pylint: disable=redefined-builtin
    FileNotFoundError = IOError

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
        providers = [ProviderBase.create(p) for 
                        p in ProviderBase.list_providers()]

    for p in providers:
        for row in p.list_available_genomes():
            yield [p.name] + list(row)

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

def list_installed_genomes(genome_dir=None):
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
    if not genome_dir:
        genome_dir = config.get("genome_dir", None)
    if not genome_dir:
        raise norns.exceptions.ConfigError("Please provide or configure a genome_dir")

    return [f for f in os.listdir(genome_dir) if 
            _is_genome_dir(genome_dir + "/" + f)]

def search(term, provider=None):
    """
    Search for a genome.

     If provider is specified, search only that specific provider, else 
     search all providers. Both the name and description are used for the 
     search. Seacrch term is case-insensitive.

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
        providers = [ProviderBase.create(p) for 
                        p in ProviderBase.list_providers()]

    for p in providers:
        for row in p.search(term):
            yield [x.encode('utf-8') for x in [p.name] + list(row)]

def install_genome(name, provider, genome_dir=None):
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
    if not genome_dir:
        genome_dir = config.get("genome_dir", None)
    if not genome_dir:
        raise norns.exceptions.ConfigError("Please provide or configure a genome_dir")
   
    genome_dir = os.path.expanduser(genome_dir)
    p = ProviderBase.create(provider)
    local_name = p.download_genome(name, genome_dir)
    generate_sizes(local_name, genome_dir)

def genome(name, genome_dir=None):
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
    if not genome_dir:
        genome_dir = config.get("genome_dir", None)
    if not genome_dir:
        raise norns.exceptions.ConfigError("Please provide or configure a genome_dir")
    
    genome_dir = os.path.expanduser(genome_dir)
    if not os.path.exists(genome_dir):
        raise FileNotFoundError(
                "genome_dir {} does not exist".format(genome_dir)
                )

    pattern = os.path.join(genome_dir, name, "*.fa")
    fnames = glob.glob(pattern)
    if len(fnames) == 0:
        raise FileNotFoundError(
                "no *.fa files found in genome_dir {}".format(
                    os.path.join(genome_dir, name)
                    )
                )
    elif len(fnames) > 1:
        fname = os.path.join(genome_dir, name, "{}.fa".format(name))
        if fname not in fnames:
            raise Exception("More than one FASTA file found, no {}.fa!".format(name))
    else:
        fname = fnames[0]

    return Fasta(fname)
