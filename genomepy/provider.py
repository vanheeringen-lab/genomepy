import sys
import requests
import re
import os
import ftplib
import urllib2
import zlib
import xmltodict

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

        if not os.path.exists(genome_dir):
            os.makedirs(genome_dir)
        
        dbname, link = self.get_genome_download_link(name, mask)
        fname = os.path.join(genome_dir, dbname, dbname + ".fa") 
        gzipped = False
        if link.endswith(".gz"):
            gzipped = True

        if not os.path.exists(os.path.join(genome_dir, dbname)):
            os.makedirs(os.path.join(genome_dir, dbname))
        response = urllib2.urlopen(link)
        
        sys.stderr.write("downloading...\n")
        with open(fname, "w") as f:
            if gzipped:
                f.write(zlib.decompress(response.read(), zlib.MAX_WBITS | 16))
            else:
                f.write(response.read())
        sys.stderr.write("done...\n")
        sys.stderr.write("name: {}\n".format(dbname))
        sys.stderr.write("fasta: {}\n".format(fname))
        

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

            response = urllib2.urlopen(self.das_url)
            d = xmltodict.parse(response.read())
            for genome in d['DASDSN']['DSN']:
                self.genomes.append([genome['SOURCE']['@id'], genome['DESCRIPTION']])
        return self.genomes

    
    def search(self, term):
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
