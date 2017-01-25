import sys
import requests
import re
import os
import ftplib
try: 
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen
import zlib
import xmltodict
import shutil
import time

from genomepy import exceptions

class ProviderBase(object):
    _providers = {}
    _name = None

    @classmethod
    def create(cls, provider):
        try:
            return cls._providers[provider]()
        except KeyError:
            raise Exception("Unknown provider")

    @classmethod
    def register_provider(cls, provider):
        def decorator(subclass):
            cls._providers[provider] = subclass
            subclass._name = provider
            return subclass
        return decorator
    
    @classmethod
    def list_providers(self):
        return self._providers.keys()

    def download_genome(self, name, genome_dir, mask="soft"):
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
        
        dbname, link = self.get_genome_download_link(name, mask)
        fname = os.path.join(genome_dir, dbname, dbname + ".fa") 
        gzipped = False
        if link.endswith(".gz"):
            gzipped = True

        if not os.path.exists(os.path.join(genome_dir, dbname)):
            os.makedirs(os.path.join(genome_dir, dbname))
        response = urlopen(link)
        
        sys.stderr.write("downloading...\n")
        with open(fname, "w") as f:
            if gzipped:
                f.write(zlib.decompress(response.read(), zlib.MAX_WBITS | 16))
            else:
                f.write(response.read())
        sys.stderr.write("done...\n")
        sys.stderr.write("name: {}\n".format(dbname))
        sys.stderr.write("fasta: {}\n".format(fname))
        
        if hasattr(self, '_post_process_download'):
            self._post_process_download(name, genome_dir)

        # Create readme with information
        readme = os.path.join(genome_dir, dbname, "README.txt")
        with open(readme, "w") as f:
            f.write("name: {}\n".format(dbname))
            f.write("original name: {}\n".format(os.path.split(link)[-1]))
            f.write("url: {}\n".format(link))
            f.write("date: {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S")))


register_provider = ProviderBase.register_provider

@register_provider('Ensembl')
class EnsemblProvider(ProviderBase):
    rest_url = "http://rest.ensemblgenomes.org/"
    
    def __init__(self):
        self.genomes = None

    def request_json(self, ext):
        r = requests.get(self.rest_url + ext, 
            headers={ "Content-Type" : "application/json"})

        if not r.ok:
            r.raise_for_status()
        
        return r.json()
    
    def list_available_genomes(self, cache=True, as_dict=False):
        ext = "/info/genomes?"
        if not cache or not self.genomes:
             
            self.genomes = self.request_json(ext)
        
        for genome in self.genomes:
            if as_dict:
                yield genome
            else:
                yield genome.get("dbname", ""), genome.get("name", "")
   
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
            if term in ",".join([str(v) for v in genome.values()]).lower():
                yield genome.get("dbname", ""), genome.get("name", "")

    def get_genome_download_link(self, name, mask="soft"):
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

        # get the genome information
        try:
            ext = "/info/genomes/" + name + "/?"
            genome_info = self.request_json(ext)
        except requests.exceptions.HTTPError as e:
            sys.stderr.write("Species not found: {}".format(e))
            raise exceptions.GenomeDownloadError(
                "Could not download genome {} from Ensembl".format(name))
        
        # parse the division
        division = genome_info["division"].lower().replace("ensembl","")
        if division == "bacteria":
            raise NotImplementedError("bacteria from ensembl not yet supported")
        
        # get version info from the dbname string
        p = re.compile(r'core_(\d+)')
        dbname = genome_info["dbname"]
        m = p.search(dbname)
        if m:
            version = m.group(1)
        species = genome_info["species"]
        
        ftp_site = "ftp.ensemblgenomes.org"
        if not division:
            ftp_site = "ftp.ensembl.org"

        if division:
            base_url = "/pub/{}/release-{}/fasta/{}/dna/"
            ftp_dir = base_url.format(division, version, species)
        else:
            base_url = "/pub/release-{}/fasta/{}/dna/"
            ftp_dir = base_url.format(version, species)
        
        ftp = ftplib.FTP(ftp_site)
        ftp.login("anonymous", "s.vanheeringen@science.ru.nl")
        fnames = ftp.nlst(ftp_dir)
        
        pattern = "dna.toplevel"
        if mask == "soft":
            pattern = "dna_sm.toplevel"
        elif mask == "hard":
            pattern = "dna_rm.toplevel"

        for fname in fnames:
            if pattern in fname:
                return dbname, "ftp://" + ftp_site + fname

@register_provider('UCSC')
class UcscProvider(ProviderBase):
    ucsc_url = "http://hgdownload.soe.ucsc.edu/goldenPath/{0}/bigZips/chromFa.tar.gz"
    alt_ucsc_url = "http://hgdownload.soe.ucsc.edu/goldenPath/{0}/bigZips/{0}.fa.gz"
    das_url = "http://genome.ucsc.edu/cgi-bin/das/dsn"
    
    def __init__(self):
        self.genomes = []
    
    def list_available_genomes(self, cache=True):
        if not cache or len(self.genomes) == 0:

            response = urlopen(self.das_url)
            d = xmltodict.parse(response.read())
            for genome in d['DASDSN']['DSN']:
                self.genomes.append([genome['SOURCE']['@id'], genome['DESCRIPTION']])
        return self.genomes
    
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
        for name,description in self.list_available_genomes():
            if term in name.lower() or term in description.lower():
                yield name,description
    
    def get_genome_download_link(self, name, mask="soft"):
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
        for genome_url in [self.ucsc_url, self.alt_ucsc_url]:
            remote = genome_url.format(name)
            ret = requests.head(remote)
            
            if ret.status_code == 200:
                return name, remote

        raise exceptions.GenomeDownloadError(
                "Could not download genome {} from UCSC".format(name))

@register_provider('NCBI')
class NCBIProvider(ProviderBase):
    assembly_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/"

    def __init__(self):
        self.genomes = None

    def _get_genomes(self):
        genomes = []
        
        names = [
                "assembly_summary_refseq.txt", 
                "assembly_summary_refseq_historical.txt",
                ]
        
        for fname in names:
            response = urlopen(self.assembly_url + "/" + fname)
            lines = response.read().decode('utf-8').splitlines()
            header = lines[1].strip("# ").split("\t")
            for line in lines[2:]:
                vals = line.strip("# ").split("\t")
                genomes.append(dict(zip(header, vals)))
        
        return genomes

    def list_available_genomes(self, cache=True, as_dict=False):
        if not cache or not self.genomes:
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

    def get_genome_download_link(self, name, mask="soft"):
        """
        Return NCBI ftp link to toplevel genome sequence

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
        
        if not self.genomes:
            self.genomes = self._get_genomes()
        
        for genome in self.genomes:
            if genome["asm_name"] == name:
                #ftp_path': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/195/GCF_000004195.3_Xenopus_tropicalis_v9.1'
                url = genome["ftp_path"]
                url += "/" + url.split("/")[-1] + "_genomic.fna.gz"
                return name, url
        raise exceptions.GenomeDownloadError("Could not download genome from NCBI")
    
    def _post_process_download(self, name, genome_dir):
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
                break
   
        # Create mapping of accessions to names
        tr = {}
        response = urlopen(url)
        for line in response.read().decode('utf-8').splitlines():
            if line.startswith("#"):
                continue
            vals = line.strip().split("\t")
            tr[vals[6]] = vals[0]
    
        # Check of the original genome fasta exists
        fa = os.path.join(genome_dir, name, "{}.fa".format(name))
        if not os.path.exists(fa):
            raise Exception("Genome fasta file not found, {}".format(fa))
        
        # Use a tmp file and replace the names
        new_fa = os.path.join(
                genome_dir, name,
                ".process.{}.fa".format(name)
                )

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
                    
                    else:
                        new.write(line)
        
        # Rename tmp file to real genome file
        shutil.move(new_fa, fa)

if __name__ == "__main__":
    p = ProviderBase.create("NCBI")
    p.download_genome("dyak_caf1", "/data/genomes")

