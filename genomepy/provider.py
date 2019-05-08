"""Genome providers."""
import sys
import requests
import re
import os
import io
import time
import gzip
import xmltodict
import shutil
import subprocess as sp
import tarfile
import time
from tempfile import mkdtemp,NamedTemporaryFile
try:
        from urllib.request import urlopen, urlretrieve, urlcleanup
except:
        from urllib import urlopen, urlretrieve, urlcleanup

from bucketcache import Bucket
from pyfaidx import Fasta
from appdirs import user_cache_dir

from genomepy import exceptions
from genomepy.utils import filter_fasta, get_localname
from genomepy.__about__ import __version__

my_cache_dir = os.path.join(user_cache_dir("genomepy"), __version__)
# Create .cache dir if it does not exist
if not os.path.exists(my_cache_dir):
    os.makedirs(my_cache_dir)

cached = Bucket(my_cache_dir, days=7)

class ProviderBase(object):
    
    """Provider base class.

    Use to get a list of available providers:
    >>> ProviderBase.list_providers()
    ['UCSC', 'NCBI', 'Ensembl']

    Create a provider:
    >>> p = ProviderBase.create("UCSC")
    >>> for name, desc in p.search("hg38"):
    ...     print(desc)
    Human Dec. 2013 (GRCh38/hg38) Genome at UCSC
    """
    
    _providers = {}
    name = None

    @classmethod
    def create(cls, name):
        """Create a provider based on the provider name.

        Parameters
        ----------
        name : str
            Name of the provider (eg. UCSC, Ensembl, ...)
        
        Returns
        -------
        provider : Provider instance
            Provider instance.
        """
        try:
            return cls._providers[name]()
        except KeyError:
            raise Exception("Unknown provider")

    @classmethod
    def register_provider(cls, provider):
        """Register method to keep list of providers."""
        def decorator(subclass):
            """Register as decorator function."""
            cls._providers[provider] = subclass
            subclass.name = provider
            return subclass
        return decorator
    
    @classmethod
    def list_providers(self):
        """List available providers.""" 
        return self._providers.keys()

    def __hash__(self):
        return hash(str(self.__class__))
    
    def tar_to_bigfile(self, fname, outfile):
        """Convert tar of multiple FASTAs to one file."""
        fnames = []
        tmpdir = mkdtemp()
        
        # Extract files to temporary directory
        with tarfile.open(fname) as tar:
            tar.extractall(path=tmpdir)
        for root, _, files in os.walk(tmpdir):
            fnames += [os.path.join(root, fname) for fname in files]
        
        # Concatenate
        with open(outfile, "w") as out:
            for infile in fnames:
                for line in open(infile):
                    out.write(line)
                os.unlink(infile)
        
        # Remove temp dir
        shutil.rmtree(tmpdir)

    def download_genome(self, name, genome_dir, localname=None, mask="soft", regex=None, invert_match=False, version=None):
        """
        Download a (gzipped) genome file to a specific directory

        Parameters
        ----------
        name : str
            Genome / species name
        
        genome_dir : str
            Directory to install genome

        mask: str , optional
            Masking, soft, hard or none (all other strings)
        """
        genome_dir = os.path.expanduser(genome_dir)
        
        if not os.path.exists(genome_dir):
            os.makedirs(genome_dir)
        
        dbname, link = self.get_genome_download_link(name, mask=mask, version=version)
        myname = dbname 
        if localname:
            myname = localname
        
        myname = myname.replace(" ", "_")

        gzipped = False
        if link.endswith(".gz"):
            gzipped = True

        if not os.path.exists(os.path.join(genome_dir, myname)):
            os.makedirs(os.path.join(genome_dir, myname))
        urlcleanup()
        response = urlopen(link)
         
        sys.stderr.write("downloading from {}...\n".format(link))
        down_dir = genome_dir
        fname = os.path.join(genome_dir, myname, myname + ".fa")
        if regex:
            down_dir = mkdtemp()
            os.mkdir(os.path.join(down_dir, myname))
            fname = os.path.join(down_dir, myname, myname + ".fa") 
        with open(fname, "wb") as f_out:
            if gzipped:
                # Supports both Python 2.7 as well as 3
                with gzip.GzipFile(fileobj=io.BytesIO(response.read())) as f_in:
                    shutil.copyfileobj(f_in, f_out)
            else:
                f_out.write(response.read())
        sys.stderr.write("done...\n")
        
        if link.endswith("tar.gz"):
            self.tar_to_bigfile(fname, fname) 
        
        if hasattr(self, '_post_process_download'):
            self._post_process_download(name, down_dir, mask)
        
        if regex:
            infa = fname
            outfa = os.path.join(genome_dir, myname, myname + ".fa") 
            filter_fasta(
                infa,
                outfa,
                regex=regex,
                v=invert_match,
                force=True
                )

            not_included = [k for k in Fasta(infa).keys() if k not in Fasta(outfa).keys()]
            shutil.rmtree(down_dir)
            fname = outfa
        
        sys.stderr.write("name: {}\n".format(dbname))
        sys.stderr.write("local name: {}\n".format(myname))
        sys.stderr.write("fasta: {}\n".format(fname))

        # Create readme with information
        readme = os.path.join(genome_dir, myname, "README.txt")
        with open(readme, "w") as f:
            f.write("name: {}\n".format(myname))
            f.write("original name: {}\n".format(dbname))
            f.write("original filename: {}\n".format(os.path.split(link)[-1]))
            f.write("url: {}\n".format(link))
            f.write("mask: {}\n".format(mask))
            f.write("date: {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S")))
            if regex:
                if invert_match:
                    f.write("regex: {} (inverted match)\n".format(regex))
                else:
                    f.write("regex: {}\n".format(regex))
                f.write("sequences that were excluded:\n")
                for seq in not_included:
                    f.write("\t{}\n".format(seq))
    
    def download_annotation(self, name, genome_dir, localname=None, version=None):
        """
        Download annotation file to to a specific directory

        Parameters
        ----------
        name : str
            Genome / species name
        
        genome_dir : str
            Directory to install annotation
        """
        raise NotImplementedError()

register_provider = ProviderBase.register_provider

@register_provider('Ensembl')
class EnsemblProvider(ProviderBase):
    
    """
    Ensembl genome provider.

    Will search both ensembl.org as well as ensemblgenomes.org.
    The bacteria division is not yet supported.
    """
    
    rest_url = "http://rest.ensembl.org/"
    
    def __init__(self):
        self.genomes = None
        self.version = None

    @cached(method=True)
    def request_json(self, ext):
        """Make a REST request and return as json."""
        if self.rest_url.endswith("/") and ext.startswith("/"):
            ext = ext[1:]
        
        r = requests.get(self.rest_url + ext, 
            headers={ "Content-Type" : "application/json"})

        if not r.ok:
            r.raise_for_status()
        
        return r.json()
    
    def safe(self, name):
        """Replace spaces with undescores.
        """
        return name.replace(" ", "_")

    def list_available_genomes(self, as_dict=False):
        """
        List all available genomes.
        
        Parameters
        ----------
        as_dict : bool, optional
            Return a dictionary of results.
        
        Yields
        ------
        genomes : dictionary or tuple
        """
        if not self.genomes:
            self.genomes = []
            divisions = self.request_json('info/divisions?')
            for division in divisions:
                genomes = self.request_json("info/genomes/division/{}?".format(division))
                self.genomes += genomes
       
        for genome in self.genomes:
            if as_dict:
                yield genome
            else:
                yield (self.safe(genome.get("assembly_name", "")), genome.get("name", ""))
   
    def search(self, term):
        """
        Search for a genome at Ensembl. 
        
        Both the name and description are used for the 
        search. Seacrch term is case-insensitive.
    
        Parameters
        ----------
        term : str
            Search term, case-insensitive.
    
        Yields
        ------
        tuple
            genome information (name/identifier and description)
        """
        term = term.lower()
        for genome in self.list_available_genomes(as_dict=True):
            if term in ",".join([str(v) for v in genome.values()]).lower():
                yield (self.safe(genome.get("assembly_name", "")), 
                        genome.get("name", ""))

    def _get_genome_info(self, name):
        """Get genome_info from json request."""
        try:
            assembly_acc = ""
            for genome in self.list_available_genomes(as_dict=True):
                if self.safe(genome.get("assembly_name", "")) == self.safe(name):
                    assembly_acc = genome.get("assembly_accession", "")
                    break
            if assembly_acc:
                ext  = "info/genomes/assembly/" + assembly_acc + "/?"
                genome_info = self.request_json(ext)
            else:
                raise exceptions.GenomeDownloadError(
                    "Could not download genome {} from Ensembl".format(name))
        except requests.exceptions.HTTPError as e:
            sys.stderr.write("Species not found: {}".format(e))
            raise exceptions.GenomeDownloadError(
                "Could not download genome {} from Ensembl".format(name))
        return genome_info

    def get_version(self, ftp_site):
        """Retrieve current version from Ensembl FTP.
        """
        response = urlopen(ftp_site + "/current_README")
        p = re.compile(r'Ensembl (Genomes|Release) (\d+)')
        m = p.search(response.read().decode())
        if m:
            version = m.group(2)
            sys.stderr.write("Using version {}\n".format(version))
            self.version = version
            return version

    def get_genome_download_link(self, name, version=None, mask="soft"):
        """
        Return Ensembl ftp link to toplevel genome sequence

        Parameters
        ----------
        name : str 
            Genome name. Current implementation will fail if exact
            name is not found.
        
        Returns
        ------
        tuple (name, link) where name is the Ensembl dbname identifier
        and link is a str with the ftp download link.
        """
        genome_info = self._get_genome_info(name)
        
        # parse the division
        division = genome_info["division"].lower().replace("ensembl","")
        if division == "bacteria":
            raise NotImplementedError("bacteria from ensembl not yet supported")
        
        ftp_site = "ftp://ftp.ensemblgenomes.org/pub"
        if division == "vertebrates":
            ftp_site = "https://ftp.ensembl.org/pub"

        version = self.version
        if not version:
            version = self.get_version(ftp_site)

        if division != 'vertebrates':
            base_url = "/{}/release-{}/fasta/{}/dna/"
            ftp_dir = base_url.format(division, version, genome_info["url_name"].lower())
            url = "{}/{}".format(ftp_site, ftp_dir)
        else:
            base_url = "/release-{}/fasta/{}/dna/"
            ftp_dir = base_url.format(version, genome_info["url_name"].lower())
            url = "{}/{}".format(ftp_site, ftp_dir)
      
        pattern = "dna.toplevel"
        if mask == "soft":
            pattern = "dna_sm.toplevel"
        elif mask == "hard":
            pattern = "dna_rm.toplevel"
        
        asm_url = "{}/{}.{}.{}.fa.gz".format(
            url,
            genome_info['url_name'].capitalize(),
            re.sub(r'\.p\d+$', '', self.safe(genome_info["assembly_name"])),
            pattern)
        
        return self.safe(genome_info["assembly_name"]), asm_url

    def download_annotation(self, name, genome_dir, localname=None, version=None):
        """
        Download gene annotation from Ensembl based on genome name.
    
        Parameters
        ----------
        name : str
            Ensembl genome name.
        genome_dir : str
            Genome directory.
        version : str , optional
            Ensembl version. By default the latest version is used.
        """
        localname = get_localname(name, localname)
        genome_info = self._get_genome_info(name)

        # parse the division
        division = genome_info["division"].lower().replace("ensembl","")
        if division == "bacteria":
            raise NotImplementedError("bacteria from ensembl not yet supported")
                
        # Get the base link depending on division
        ftp_site = "https://ftp.ensemblgenomes.org/pub/{}".format(division)
        if division == 'vertebrates':
            ftp_site = "https://ftp.ensembl.org/pub"

        if not version:
            version = self.get_version(ftp_site)

        # Get the GTF URL
        base_url = ftp_site + "/release-{}/gtf/{}/{}.{}.{}.gtf.gz"
        safe_name = name.replace(" ", "_")
        safe_name = re.sub('\.p\d+$', '', safe_name)

        ftp_link = base_url.format(version, genome_info["url_name"].lower(), genome_info["url_name"].capitalize(), safe_name, version)

        out_dir = os.path.join(genome_dir, localname)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        # Download the file
        try:
            response = urlopen(ftp_link)
            gtf_file = out_dir + "/"+ localname + ".annotation.gtf.gz"
            with open(gtf_file, "wb") as f:
                f.write(response.read())
        
            bed_file = gtf_file.replace("gtf.gz", "bed")
            cmd = "gtfToGenePred {0} /dev/stdout | genePredToBed /dev/stdin {1} && gzip {1}"
            ret = sp.check_call(cmd.format(gtf_file, bed_file), shell=True)
            readme = os.path.join(genome_dir, localname, "README.txt")
            with open(readme, "a") as f:
                f.write("annotation url: {}\n".format(ftp_link))
        except Exception:
            sys.stderr.write("\nCould not download {}\n".format(ftp_link))
            raise

        return out_dir

@register_provider('UCSC')
class UcscProvider(ProviderBase):
     
    """
    UCSC genome provider.
    
    The UCSC DAS server is used to search and list genomes. 
    """
    
    ucsc_url = "http://hgdownload.soe.ucsc.edu/goldenPath/{0}/bigZips/chromFa.tar.gz"
    ucsc_url_masked = "http://hgdownload.soe.ucsc.edu/goldenPath/{0}/bigZips/chromFaMasked.tar.gz"
    alt_ucsc_url = "http://hgdownload.soe.ucsc.edu/goldenPath/{0}/bigZips/{0}.fa.gz"
    alt_ucsc_url_masked = "http://hgdownload.soe.ucsc.edu/goldenPath/{0}/bigZips/{0}.fa.masked.gz"
    das_url = "http://genome.ucsc.edu/cgi-bin/das/dsn"
    
    def __init__(self):
        self.genomes = []
   
    def list_available_genomes(self):
        """
        List all available genomes.
        
        Returns
        -------
        genomes : list
        """
        self.genomes = self._get_genomes()
        
        return self.genomes
    
    @cached(method=True) 
    def _get_genomes(self):
        response = urlopen(self.das_url)
        d = xmltodict.parse(response.read())
        genomes = []
        for genome in d['DASDSN']['DSN']:
                genomes.append([
                    genome['SOURCE']['@id'], genome['DESCRIPTION']
                    ])
        return genomes


    def search(self, term):
        """
        Search for a genome at UCSC. 
        
        Both the name and description are used for the 
        search. Search term is case-insensitive.
    
        Parameters
        ----------
        term : str
            Search term, case-insensitive.
    
        Yields
        ------
        tuple
            genome information (name/identifier and description)
        """
        term = term.lower()
        for name,description in self.list_available_genomes():
            if term in name.lower() or term in description.lower():
                yield name,description
    
    def get_genome_download_link(self, name, mask="soft", version=None):
        """
        Return UCSC http link to genome sequence

        Parameters
        ----------
        name : str 
            Genome build name. Current implementation will fail if exact
            name is not found.
        
        Returns
        ------
        tuple (name, link) where name is the genome build identifier
        and link is a str with the http download link.
        """
        urls = [self.ucsc_url, self.alt_ucsc_url]
        if mask == "hard":
            urls = [self.ucsc_url_masked, self.alt_ucsc_url_masked]

        for genome_url in urls:
            remote = genome_url.format(name)
            ret = requests.head(remote)
            
            if ret.status_code == 200:
                return name, remote

        raise exceptions.GenomeDownloadError(
                "Could not download genome {} from UCSC".format(name))


    def download_annotation(self, name, genome_dir, localname=None, version=None):
        """
        Download gene annotation from UCSC based on genomebuild.
    
        Will check UCSC, Ensembl and RefSeq annotation.
    
        Parameters
        ----------
        genomebuild : str
            UCSC genome name.
        genome_dir : str
            Genome directory.
        """
        localname = get_localname(name, localname)

        UCSC_GENE_URL = "http://hgdownload.cse.ucsc.edu/goldenPath/{}/database/"
        ANNOS = ["knownGene.txt.gz", "ensGene.txt.gz", "refGene.txt.gz"]
        pred = "genePredToBed"
    
        tmp = NamedTemporaryFile(delete=False, suffix=".gz")
    
        anno = []
        f = urlopen(UCSC_GENE_URL.format(name))
        p = re.compile(r'\w+.Gene.txt.gz')
        for line in f.readlines():
            m = p.search(line.decode())
            if m:
                anno.append(m.group(0))
        sys.stderr.write("Retrieving gene annotation for {}\n".format(name))
        url = ""
        for a in ANNOS:
            if a in anno:
                url = UCSC_GENE_URL.format(name) + a
                break
        if url:
            sys.stderr.write("Using {}\n".format(url))
            urlretrieve(
                    url,
                    tmp.name
                    )
             
            with gzip.open(tmp.name) as f:
                cols = f.readline().decode(errors='ignore').split("\t")
    
            start_col = 1
            for i,col in enumerate(cols):
                if col == "+" or col == "-":
                    start_col = i - 1
                    break
            end_col = start_col + 10
           
            localname = localname.replace(" ", "_")
            path = os.path.join(genome_dir, localname)
            if not os.path.exists(path):
                os.mkdir(path)
            bed_file = os.path.join(genome_dir, localname, localname + ".annotation.bed")
            cmd = "zcat {} | cut -f{}-{} | {} /dev/stdin {} && gzip {}"
            sp.call(cmd.format(
                tmp.name, start_col, end_col, pred, bed_file, bed_file), 
                shell=True)
            
            gtf_file = bed_file.replace(".bed", ".gtf")
            cmd = "bedToGenePred {0}.gz /dev/stdout | genePredToGtf file /dev/stdin /dev/stdout -utr -honorCdsStat | sed 's/.dev.stdin/UCSC/' > {1} && gzip {1}"
            ret = sp.check_call(cmd.format(bed_file, gtf_file), shell=True)
        
            readme = os.path.join(genome_dir, localname, "README.txt")
            with open(readme, "a") as f:
                f.write("annotation url: {}\n".format(url))
        else:
            sys.stderr.write("No annotation found!")

@register_provider('NCBI')
class NCBIProvider(ProviderBase):
        
    """
    NCBI genome provider.
    
    Useas the assembly reports page to search and list genomes. 
    """
    
    assembly_url = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/"

    def __init__(self):
        self.genomes = None

    @cached(method=True)
    def _get_genomes(self):
        """Parse genomes from assembly summary txt files."""
        genomes = []
        
        names = [
                "assembly_summary_refseq.txt", 
                "assembly_summary_genbank.txt", 
                "assembly_summary_refseq_historical.txt",
                ]
        
        sys.stderr.write("Downloading assembly summaries from NCBI, this will take a while...\n")
        seen = {}
        for fname in names:
            urlcleanup() 
            response = urlopen(self.assembly_url + "/" + fname)
            lines = response.read().decode('utf-8').splitlines()
            header = lines[1].strip("# ").split("\t")
            for line in lines[2:]:
                vals = line.strip("# ").split("\t")
                # Don't repeat samples with the same BioSample ID
                if vals[2] not in seen:   # BioSample ID
                    genomes.append(dict(zip(header, vals)))
                    seen[vals[2]] = 1
        
        return genomes

    def list_available_genomes(self, as_dict=False):
        """
        List all available genomes.
        
        Parameters
        ----------
        as_dict : bool, optional
            Return a dictionary of results.
        
        Yields
        ------
        genomes : dictionary or tuple
        """
        if not self.genomes:
            self.genomes = self._get_genomes()

        for genome in self.genomes:
            if as_dict:
                yield genome
            else:
                yield (
                        genome.get("asm_name", ""), 
                        "; ".join((
                            genome.get("organism_name", ""),
                            genome.get("submitter", ""),
                        ))
                    )
   
    def search(self, term):
        """
        Search for term in genome names and descriptions. 

        The search is case-insensitive.

        Parameters
        ----------
        term : str
            search term
        
        Yields
        ------
        tuples with two items, name and description
        """
        term = term.lower()

        for genome in self.list_available_genomes(as_dict=True):
            term_str = ";".join([repr(x) for x in genome.values()])

            if term in term_str.lower():
                yield (
                        genome.get("asm_name", ""), 
                        "; ".join((
                            genome.get("organism_name", ""),
                            genome.get("submitter", ""),
                        ))
                    )

    def get_genome_download_link(self, name, mask="soft", version=None):
        """
        Return NCBI ftp link to top-level genome sequence

        Parameters
        ----------
        name : str 
            Genome name. Current implementation will fail if exact
            name is not found.
        
        Returns
        ------
        tuple (name, link) where name is the NCBI asm_name identifier
        and link is a str with the ftp download link.
        """
        if mask == "hard":
            sys.stderr.write("ignoring mask parameter for NCBI at download.\n")

        if not self.genomes:
            self.genomes = self._get_genomes()
        
        for genome in self.genomes:
            if genome["asm_name"] == name:
                #ftp_path': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/195/GCF_000004195.3_Xenopus_tropicalis_v9.1'
                url = genome["ftp_path"]
                url = url.replace("ftp://", "https://")
                url += "/" + url.split("/")[-1] + "_genomic.fna.gz"
                return name, url
        raise exceptions.GenomeDownloadError("Could not download genome from NCBI")
    
    def _post_process_download(self, name, genome_dir, mask="soft"):
        """
        Replace accessions with sequence names in fasta file.
        
        Parameters
        ----------
        name : str
            NCBI genome name

        genome_dir : str
            Genome directory
        """
        # Get the FTP url for this specific genome and download
        # the assembly report
        for genome in self.genomes:
            if genome["asm_name"] == name:
                url = genome["ftp_path"]
                url += "/" + url.split("/")[-1] + "_assembly_report.txt"
                url = url.replace("ftp://", "https://")
                break
   
        # Create mapping of accessions to names
        tr = {}
        urlcleanup()
        response = urlopen(url)
        for line in response.read().decode('utf-8').splitlines():
            if line.startswith("#"):
                continue
            vals = line.strip().split("\t")
            tr[vals[6]] = vals[0]
    
        name = name.replace(" ", "_")
        # Check of the original genome fasta exists
        fa = os.path.join(genome_dir, name, "{}.fa".format(name))
        if not os.path.exists(fa):
            raise Exception("Genome fasta file not found, {}".format(fa))
        
        # Use a tmp file and replace the names
        new_fa = os.path.join(
                genome_dir, name,
                ".process.{}.fa".format(name)
                )
        if mask == "hard":
            sys.stderr.write("masking lower-case.\n")

        with open(fa) as old:
            with open(new_fa, "w") as new:
                for line in old:
                    if line.startswith(">"):
                        desc = line.strip()[1:]
                        name = desc.split(" ")[0]
                        new.write(">{} {}\n".format(
                            tr.get(name, name),
                            desc
                            ))
                    elif mask == "hard":
                        new.write(re.sub('[actg]', 'N', line))
                    else:
                        new.write(line)
        
        # Rename tmp file to real genome file
        shutil.move(new_fa, fa)

    def download_annotation(self, name, genome_dir, localname=None, version=None):
        """
        Download annotation file to to a specific directory

        Parameters
        ----------
        name : str
            Genome / species name
        
        genome_dir : str
            Directory to install annotation
        """
        sys.stderr.write("Downloading annotation...\n")
        if not self.genomes:
            self.genomes = self._get_genomes()
        
        for genome in self.genomes:
            if genome["asm_name"] == name:
                #ftp_path': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/195/GCF_000004195.3_Xenopus_tropicalis_v9.1'
                url = genome["ftp_path"]
                url += "/" + url.split("/")[-1] + "_genomic.gff.gz"
  
        localname = get_localname(name, localname)
        out_dir = os.path.join(genome_dir, localname)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        
        # Download the file
        try:
            response = urlopen(url)
            gff_file = out_dir + "/"+ localname + ".annotation.gff.gz"
            with open(gff_file, "wb") as f:
                f.write(response.read())
        except Exception:
            sys.stderr.write("WARNING: Could not download annotation from NCBI, skipping.\n")
            sys.stderr.write("URL: {}\n".format(url))

            sys.stderr.write("If you think the annotation should be there, please file a bug report at:\n")
            sys.stderr.write("https://github.com/simonvh/genomepy/issues\n")
            return
        
        cmd = "gff3ToGenePred {0} /dev/stdout | wc -l"
        out = sp.check_output(cmd.format(gff_file), shell=True)
        if out.strip() == b"0":
            sys.stderr.write("WARNING: annotation from NCBI contains no genes, skipping.\n")
        else:
            # Convert to BED file
            bed_file = gff_file.replace("gff.gz", "bed")
            cmd = "gff3ToGenePred -rnaNameAttr=gene {0} /dev/stdout | genePredToBed /dev/stdin {1} && gzip {1}"
            ret = sp.check_call(cmd.format(gff_file, bed_file), shell=True) 

            # Convert to GTF file
            gtf_file = gff_file.replace("gff.gz", "gtf")
            cmd = "gff3ToGenePred -geneNameAttr=gene {0} /dev/stdout | genePredToGtf file /dev/stdin {1} && gzip {1}"
            ret = sp.check_call(cmd.format(gff_file, gtf_file), shell=True)
        
        readme = os.path.join(genome_dir, name, "README.txt")
        with open(readme, "a") as f:
            f.write("annotation url: {}\n".format(url))

        return out_dir



