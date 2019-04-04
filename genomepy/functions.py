"""Module-level functions."""
import bisect
import os
import glob
import random
import re
import shutil
from tempfile import mkdtemp

from appdirs import user_config_dir
from pyfaidx import Fasta,Sequence
from genomepy.provider import ProviderBase
from genomepy.plugin import get_active_plugins, init_plugins, activate
from genomepy.utils import get_localname
import norns

config = norns.config("genomepy", default="cfg/default.yaml")

# Python 2
try:
    FileNotFoundError
except NameError:
    # pylint: disable=redefined-builtin
    FileNotFoundError = IOError

def manage_config(cmd, *args):
    """Manage genomepy config file."""
    if cmd == "file":
        print(config.config_file)
    elif cmd == "show":
        with open(config.config_file) as f:
            print(f.read())
    elif cmd == "generate":
        fname = os.path.join(
                user_config_dir("genomepy"), "{}.yaml".format("genomepy")
            )
        
        if not os.path.exists(user_config_dir("genomepy")):
            os.makedirs(user_config_dir("genomepy"))
        
        with open(fname, "w") as fout:
            with open(config.config_file) as fin:
                fout.write(fin.read())
        print("Created config file {}".format(fname))

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
            yield [x.encode('latin-1') for x in [p.name] + list(row)]

def install_genome(name, provider, version=None, genome_dir=None, localname=None, mask="soft", regex=None, invert_match=False, annotation=False):
    """
    Install a genome.

    Parameters
    ----------
    name : str
        Genome name

    provider : str
        Provider name

    version : str
        Version (only for Ensembl)

    genome_dir : str , optional
        Where to store the fasta files
    
    localname : str , optional
        Custom name for this genome.

    mask : str , optional
        Default is 'soft', specify 'hard' for hard masking.

    regex : str , optional
        Regular expression to select specific chromosome / scaffold names.

    invert_match : bool , optional
        Set to True to select all chromosomes that don't match the regex.

    annotation : bool , optional
        If set to True, download gene annotation in BED and GTF format.
    """
    if not genome_dir:
        genome_dir = config.get("genome_dir", None)
    if not genome_dir:
        raise norns.exceptions.ConfigError("Please provide or configure a genome_dir")
   
    genome_dir = os.path.expanduser(genome_dir)
    localname = get_localname(name, localname)
    
    # Download genome from provider
    p = ProviderBase.create(provider)
    p.download_genome(
        name, 
        genome_dir, 
        version=version,
        mask=mask, 
        localname=localname, 
        regex=regex, 
        invert_match=invert_match)

    if annotation:
        # Download annotation from provider
        p.download_annotation(name, genome_dir, localname=localname, version=version)

    g = Genome(localname, genome_dir=genome_dir)
    for plugin in get_active_plugins():
        plugin.after_genome_download(g)

    generate_env()

def get_track_type(track):
    region_p = re.compile(r'^(.+):(\d+)-(\d+)$')
    if type(track) == type([]):
        if region_p.search(track[0]):
            return "interval"

    with open(track) as fin:
        line = fin.readline().strip()
    if region_p.search(line):
        return "interval"
    return "bed"

def _weighted_selection(l, n):
    """
        Selects  n random elements from a list of (weight, item) tuples.
        Based on code snippet by Nick Johnson
    """
    cuml = []
    items = []
    total_weight = 0.0
    for weight, item in l:
        total_weight += weight
        cuml.append(total_weight)
        items.append(item)

    return [items[bisect.bisect(cuml, random.random()*total_weight)] for _ in range(n)]

def generate_exports():
	"""Print export commands for setting environment variables.
	"""
	env = []
	for name in list_installed_genomes():
	    try:
	        g = Genome(name)
	        env_name = re.sub(r'[^\w]+', "_", name).upper()
	        env.append("export {}={}".format(env_name, g.filename))
	    except:
	        pass
	return env

def generate_env(fname=None):
    """Generate file with exports. 

    By default this is in .config/genomepy/exports.txt.

    Parameters
    ----------
    fname: strs, optional
        Name of the output file.
    """
    config_dir = user_config_dir("genomepy")
    if os.path.exists(config_dir):
        fname = os.path.join(config_dir, "exports.txt")
        with open(fname, "w") as fout:
            for env in generate_exports():
                fout.write("{}\n".format(env))

class Genome(Fasta):
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
    def __init__(self, name, genome_dir=None):
        
        try:
            super(Genome, self).__init__(name)
            self.name = os.path.basename(name)
        except:
            if os.path.isdir(name) and len(glob.glob(os.path.join(name, "*.fa"))) == 1 and genome_dir is not None:
                fname = glob.glob(os.path.join(name, "*.fa"))[0]
                name = os.path.basename(fname)  
            else:
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
    
            super(Genome, self).__init__(fname)
        
            self.name = name
        
        self._gap_sizes = None
        self.props = {}

        for plugin in get_active_plugins():
            self.props[plugin.name()] = plugin.get_properties(self) 

    def _bed_to_seqs(self, track, stranded=False, extend_up=0, extend_down=0):
        BUFSIZE = 10000
        with open(track) as fin:
            lines = fin.readlines(BUFSIZE)
            while lines:
                for line in lines:
                    if line.startswith("#") or line.startswith("track"):
                        continue
    
                    vals = line.strip().split("\t")
                    try:
                        start, end = int(vals[1]), int(vals[2])
                    except ValueError:
                        raise
    
                    rc = False
                    if stranded:
                        try:
                            rc = vals[5] == "-"
                        except IndexError:
                            pass
    
                    starts = [start]
                    ends = [end]
                    
                    chrom = vals[0]
                    
                    # BED12
                    if len(vals) == 12:
                        starts = [int(x) for x in vals[11].split(",")[:-1]]
                        sizes = [int(x) for x in vals[10].split(",")[:-1]]
                        starts = [start + x  for x in starts]
                        ends = [start + size  for start,size in zip(starts, sizes)]
                    name = "{}:{}-{}".format(chrom, start, end)    
                    try:
                        name = " ".join((name, vals[3]))
                    except: 
                        pass

                    starts = [start + 1 for start in starts]
                    
                    # extend
                    if extend_up:
                        if rc:
                            ends[-1] += extend_up
                        else:
                            starts[0] -= extend_up
                    if extend_down:
                        if rc:
                            starts[0] -= extend_down
                        else:
                            ends[-1] += extend_down
    
                    intervals = zip(starts, ends)
                    seq = self.get_spliced_seq(chrom, intervals, rc)
                    yield Sequence(name, seq.seq)
                
                lines = fin.readlines(BUFSIZE)

    def _region_to_seqs(self, track, extend_up=0, extend_down=0):
        BUFSIZE = 10000
        if type(track) == type([]):
            for name in track:
                    chrom, coords = name.split(":")
                    start, end = [int(c) for c in coords.split("-")]
                    start += 1
                    start -= extend_up
                    end += extend_down
                    seq = self.get_seq(chrom, start, end)
                    yield Sequence(name, seq.seq)
        else:
            with open(track) as fin:
                lines = fin.readlines(BUFSIZE)
                while lines:
                    for line in lines:
                        name = line.strip()
                        chrom, coords = name.split(":")
                        start, end = [int(c) for c in coords.split("-")]
                        start += 1
                        start -= extend_up
                        end += extend_down
                        seq = self.get_seq(chrom, start, end)
                        yield Sequence(name, seq.seq)
    
                    lines = fin.readlines(BUFSIZE)
    
    
    
    def track2fasta(self, track, fastafile=None, stranded=False, extend_up=0, extend_down=0):
        track_type = get_track_type(track)
        if track_type == "interval":
            seqqer = self._region_to_seqs(track, extend_up=extend_up, extend_down=extend_down)
        else:
            seqqer = self._bed_to_seqs(
                    track, 
                    stranded=stranded, 
                    extend_up=extend_up, 
                    extend_down=extend_down)
        
        if fastafile:
            with open(fastafile, "w") as fout:
                for seq in seqqer:
                    fout.write("{}\n".format(seq.__repr__()))
        else:
            return [seq for seq in seqqer]

    def gap_sizes(self):
        """Return gap sizes per chromosome.
        
        Returns
        -------
        gap_sizes : dict
            a dictionary with chromosomes as key and the total number of
            Ns as values
        """ 
        if not self._gap_sizes:
            gap_file = self.props["gaps"]["gaps"]
            self._gap_sizes = {}
            with open(gap_file) as f:
                for line in f:
                    chrom, start, end = line.strip().split("\t")
                    start, end = int(start), int(end)
                    self._gap_sizes[chrom] = self._gap_sizes.get(chrom, 0) + end - start
        return self._gap_sizes
    
    def get_random_sequences(self, n=10, length=200, chroms=None, max_n=0.1):
        """Return random genomic sequences.

        Parameters
        ----------
        n : int , optional
            Number of sequences to return.

        length : int , optional
            Length of sequences to return.

        chroms : list , optional
            Return sequences only from these chromosomes.

        max_n : float , optional
            Maximum fraction of Ns.

        Returns
        -------
        coords : list
            List with [chrom, start, end] genomic coordinates.
        """ 
        retries = 100
        cutoff = length * max_n
        if not chroms:
            chroms = self.keys()
   
        try:
            gap_sizes = self.gap_sizes()
        except:
            gap_sizes = {}
        sizes = dict([(chrom, len(self[chrom]) - gap_sizes.get(chrom, 0)) for chrom in chroms])
   
        l = [(sizes[x], x) for x in chroms if 
                sizes[x] / len(self[x]) > 0.1 and sizes[x] > 10 * length]
        chroms = _weighted_selection(l, n)
        coords = []
        
        count = {}
        for chrom in chroms:
            if chrom in count:
                count[chrom] += 1
            else:
                count[chrom] = 1

        for chrom in chroms:
            for i in range(retries):
                start = int(random.random() * (sizes[chrom] - length))
                end = start + length
                count_n = self[chrom][start:end].seq.upper().count("N")
                if count_n <= cutoff:
                    break
            if count_n > cutoff:
                raise ValueError("Failed to find suitable non-N sequence for {}".format(chrom))
            
            coords.append([chrom, start, end])

        return coords

def manage_plugins(command, plugin_names=None):
    """Enable or disable plugins.
    """ 
    if plugin_names is None:
        plugin_names = []
    active_plugins = config.get("plugin", [])
    plugins = init_plugins()
    if command == "enable":
        for name in plugin_names:
            if name not in plugins:
                raise ValueError("Unknown plugin: {}".format(name))
            if name not in active_plugins:
                active_plugins.append(name)
    elif command == "disable":
        for name in plugin_names:
            if name in active_plugins:
                active_plugins.remove(name)
    elif command == "list":
        print("{:20}{}".format("plugin", "enabled"))
        for plugin in sorted(plugins):
            print("{:20}{}".format(plugin, {False:"", True:"*"}[plugin in active_plugins]))
    else:
        raise ValueError("Invalid plugin command")
    config["plugin"] = active_plugins
    config.save()

    if command in ["enable", "disable"]:
        print("Enabled plugins: {}".format(", ".join(sorted(active_plugins))))
