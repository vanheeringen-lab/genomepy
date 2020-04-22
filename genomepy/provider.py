"""Genome providers."""
import sys
import requests
import re
import os
import norns
import time
import shutil
import subprocess as sp

from psutil import virtual_memory
from tempfile import TemporaryDirectory
from urllib.request import urlopen, urlretrieve, urlcleanup
from bucketcache import Bucket
from pyfaidx import Fasta
from appdirs import user_cache_dir

from genomepy.exceptions import GenomeDownloadError
from genomepy.utils import (
    read_readme,
    write_readme,
    filter_fasta,
    get_localname,
    tar_to_bigfile,
    get_file_info,
    read_url,
    safe,
    check_url,
    is_number,
    mkdir_p,
)
from genomepy.__about__ import __version__

# store the output of slow commands (marked with @cached) for fast reuse
my_cache_dir = os.path.join(user_cache_dir("genomepy"), __version__)
if not os.path.exists(my_cache_dir):
    os.makedirs(my_cache_dir)
cached = Bucket(my_cache_dir, days=7)

config = norns.config("genomepy", default="cfg/default.yaml")


class ProviderBase(object):
    """
    Provider base class.

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
    # class variables set by child classes:
    name = None
    genomes = {}
    accession_fields = []
    taxid_fields = []
    description_fields = []

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
            return cls._providers[name.lower()]()
        except KeyError:
            raise ValueError("Unknown provider")

    @classmethod
    def register_provider(cls, provider):
        """Register method to keep list of providers."""

        def decorator(subclass):
            """Register as decorator function."""
            cls._providers[provider.lower()] = subclass
            subclass.name = provider.lower()
            return subclass

        return decorator

    @classmethod
    def list_providers(cls):
        """List available providers."""
        return cls._providers.keys()

    def __hash__(self):
        return hash(str(self.__class__))

    def list_install_options(self, name=None):
        """List provider specific install options"""
        if name is None:
            return {}
        elif name.lower() not in self._providers:
            raise ValueError(f"Unknown provider: {name}")
        else:
            provider = self._providers[name.lower()]
            return provider.list_install_options(self)

    def _genome_info_tuple(self, name):
        """tuple with assembly metadata"""
        raise NotImplementedError()

    def list_available_genomes(self):
        """
        List all available genomes.

        Yields
        ------
        genomes : list of tuples
            tuples with assembly name, accession, scientific_name, taxonomy id and description
        """
        for genome in self.genomes:
            yield self._genome_info_tuple(genome)

    def check_name(self, name):
        """check if genome name can be found for provider"""
        if self.name == "URL":
            return

        if not safe(name) in self.genomes:
            raise GenomeDownloadError(
                f"Could not download genome {name} from {self.name}.\n\n"
                "Check for typos or try\n"
                f"  genomepy search {name} -p {self.name}\n"
            )

    def get_genome_download_link(self, name, mask="soft", **kwargs):
        raise NotImplementedError()

    def genome_taxid(self, genome):
        """Return the taxonomy_id for a genome.

        Parameters
        ----------
        genome : dict
            provider metadata dict of a genome.

        Returns
        ------
        Taxonomy id : int
        """
        for field in self.taxid_fields:
            tid = genome.get(field)
            if is_number(tid):
                return int(tid)
        return 0

    def assembly_accession(self, genome):
        """Return the assembly accession (GCA_*) for a genome.

        Parameters
        ----------
        genome : dict
            provider metadata dict of a genome.

        Returns
        ------
        str
            Assembly accession.
        """
        for key in self.accession_fields:
            accession = str(genome.get(key))
            if accession.startswith("GCA"):
                return accession
        return "na"

    def download_genome(
        self,
        name,
        genomes_dir,
        localname=None,
        mask="soft",
        regex=None,
        invert_match=False,
        bgzip=None,
        **kwargs,
    ):
        """
        Download a (gzipped) genome file to a specific directory

        Parameters
        ----------
        name : str
            Genome / species name

        genomes_dir : str
            Directory to install genome

        localname : str , optional
            Custom name for your genome

        mask: str , optional
            Masking, soft, hard or none (all other strings)

        regex : str , optional
            Regular expression to select specific chromosome / scaffold names.

        invert_match : bool , optional
            Set to True to select all chromosomes that don't match the regex.

        bgzip : bool , optional
            If set to True the genome FASTA file will be compressed using bgzip.
            If not specified, the setting from the configuration file will be used.
        """
        self.check_name(name)

        link = self.get_genome_download_link(name, mask=mask, **kwargs)

        original_name = name
        name = safe(name)
        localname = get_localname(name, localname)

        genomes_dir = os.path.expanduser(genomes_dir)
        out_dir = os.path.join(genomes_dir, localname)
        if not os.path.exists(out_dir):
            mkdir_p(out_dir)

        sys.stderr.write(f"Downloading genome from {link}...\n")

        # download to tmp dir. Move genome on completion.
        # tmp dir is in genome_dir to prevent moving the genome between disks
        with TemporaryDirectory(dir=out_dir) as tmp_dir:
            fname = os.path.join(tmp_dir, f"{localname}.fa")

            # actual download
            urlcleanup()
            with urlopen(link) as response:
                # check available memory vs file size.
                available_memory = int(virtual_memory().available)
                file_size = int(response.info()["Content-Length"])
                # download file in chunks if >75% of memory would be used
                cutoff = int(available_memory * 0.75)
                chunk_size = None if file_size < cutoff else cutoff
                with open(fname, "wb") as f_out:
                    shutil.copyfileobj(response, f_out, chunk_size)
            sys.stderr.write(
                "Genome download successful, starting post processing...\n"
            )

            # unzip genome
            if link.endswith(".tar.gz"):
                tar_to_bigfile(fname, fname)
            elif link.endswith(".gz"):
                os.rename(fname, fname + ".gz")
                ret = sp.check_call(["gunzip", "-f", fname])
                if ret != 0:
                    raise Exception(f"Error gunzipping genome {fname}")

            # process genome (e.g. masking)
            if hasattr(self, "_post_process_download"):
                self._post_process_download(
                    name=name, localname=localname, out_dir=tmp_dir, mask=mask
                )

            if regex:
                os.rename(fname, fname + "_to_regex")
                infa = fname + "_to_regex"
                outfa = fname
                filter_fasta(infa, outfa, regex=regex, v=invert_match, force=True)

                not_included = [
                    k for k in Fasta(infa).keys() if k not in Fasta(outfa).keys()
                ]

            # bgzip genome if requested
            if bgzip or config.get("bgzip"):
                ret = sp.check_call(["bgzip", "-f", fname])
                if ret != 0:
                    raise Exception(f"Error bgzipping {name}. Is tabix installed?")
                fname += ".gz"

            # transfer the genome from the tmpdir to the genome_dir
            src = fname
            dst = os.path.join(genomes_dir, localname, os.path.basename(fname))
            shutil.move(src, dst)

        sys.stderr.write("\n")
        sys.stderr.write("name: {}\n".format(name))
        sys.stderr.write("local name: {}\n".format(localname))
        sys.stderr.write("fasta: {}\n".format(dst))

        # Create readme with information
        readme = os.path.join(genomes_dir, localname, "README.txt")
        metadata = {
            "name": localname,
            "provider": self.name,
            "original name": original_name,
            "original filename": os.path.split(link)[-1],
            "assembly_accession": self.assembly_accession(self.genomes.get(name)),
            "tax_id": self.genome_taxid(self.genomes.get(name)),
            "mask": mask,
            "genome url": link,
            "annotation url": "na",
            "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        lines = []
        if regex:
            regex_line = f"regex: {regex}"
            if invert_match:
                regex_line += " (inverted match)"
            lines += ["", regex_line, "sequences that were excluded:"]
            for seq in not_included:
                lines.append(f"\t{seq}")
        write_readme(readme, metadata, lines)

    def get_annotation_download_link(self, name, **kwargs):
        raise NotImplementedError()

    @staticmethod
    def download_and_generate_annotation(genomes_dir, annot_url, localname):
        """download annotation file, convert to intermediate file and generate output files"""

        # create output directory if missing
        out_dir = os.path.join(genomes_dir, localname)
        if not os.path.exists(out_dir):
            mkdir_p(out_dir)

        # download to tmp dir. Move files on completion.
        with TemporaryDirectory(dir=out_dir) as tmpdir:
            ext, gz = get_file_info(annot_url)
            annot_file = os.path.join(tmpdir, localname + ".annotation" + ext)
            urlretrieve(annot_url, annot_file)

            # unzip input file (if needed)
            if gz:
                cmd = "mv {0} {1} && gunzip -f {1}"
                sp.check_call(cmd.format(annot_file, annot_file + ".gz"), shell=True)

            # generate intermediate file (GenePred)
            pred_file = annot_file.replace(ext, ".gp")
            if "bed" in ext:
                cmd = "bedToGenePred {0} {1}"
            elif "gff" in ext:
                cmd = "gff3ToGenePred -geneNameAttr=gene {0} {1}"
            elif "gtf" in ext:
                cmd = "gtfToGenePred {0} {1}"
            elif "txt" in ext:
                # UCSC annotations only
                with open(annot_file) as f:
                    cols = f.readline().split("\t")

                start_col = 1
                for i, col in enumerate(cols):
                    if col == "+" or col == "-":
                        start_col = i - 1
                        break
                end_col = start_col + 10
                cmd = f"cat {{0}} | cut -f{start_col}-{end_col} > {{1}}"
            else:
                raise TypeError(f"file type extension {ext} not recognized!")

            sp.check_call(cmd.format(annot_file, pred_file), shell=True)

            # generate gzipped gtf file (if required)
            gtf_file = annot_file.replace(ext, ".gtf")
            if "gtf" not in ext:
                cmd = "genePredToGtf file {0} {1} && gzip -f {1}"
                sp.check_call(cmd.format(pred_file, gtf_file), shell=True)

            # generate gzipped bed file (if required)
            bed_file = annot_file.replace(ext, ".bed")
            if "bed" not in ext:
                cmd = "genePredToBed {0} {1} && gzip -f {1}"
                sp.check_call(cmd.format(pred_file, bed_file), shell=True)

            # if input file was gtf/bed, gzip it
            if ext in [".gtf", ".bed"]:
                cmd = "gzip -f {}"
                sp.check_call(cmd.format(annot_file), shell=True)

            # transfer the files from the tmpdir to the genome_dir
            for f in [gtf_file + ".gz", bed_file + ".gz"]:
                src = f
                dst = os.path.join(out_dir, os.path.basename(f))
                shutil.move(src, dst)

    def attempt_and_report(self, name, localname, link, genomes_dir):
        if not link:
            sys.stderr.write(
                f"Could not download genome annotation for {name} from {self.name}.\n"
            )
            return

        sys.stderr.write(f"\nDownloading annotation from {link}...\n")
        try:
            self.download_and_generate_annotation(genomes_dir, link, localname)
        except Exception:
            raise GenomeDownloadError(
                f"\nCould not download annotation for {name} from {self.name}\n"
                "If you think the annotation should be there, please file a bug report at:\n"
                "https://github.com/vanheeringen-lab/genomepy/issues\n"
            )

        # TODO sanity check for genes
        sys.stderr.write(f"Annotation download successful\n")

        # Update readme annotation URL, or make a new
        readme = os.path.join(genomes_dir, localname, "README.txt")
        metadata, lines = read_readme(readme)
        metadata["annotation url"] = link
        write_readme(readme, metadata, lines)

    def download_annotation(self, name, genomes_dir, localname=None, **kwargs):
        """
        Download annotation file to to a specific directory

        Parameters
        ----------
        name : str
            Genome / species name

        genomes_dir : str
            Directory to install annotation

        localname : str , optional
            Custom name for your genome
        """
        self.check_name(name)

        link = self.get_annotation_download_link(name, **kwargs)

        self.attempt_and_report(name, localname, link, genomes_dir)

    def _search_taxids(self, genome, term):
        """check if search term corresponds to the provider's taxonomy field(s)"""
        for field in self.taxid_fields:
            if term == str(genome[field]):
                return True

    def _search_descriptions(self, genome, term):
        """check if search term corresponds to the provider's description field(s)"""
        for field in self.description_fields:
            if term in safe(genome[field].lower()):
                return True

    def search(self, term):
        """
        Search for term in genome names, descriptions and taxonomy ID.

        The search is case-insensitive.

        Parameters
        ----------
        term : str
            Search term, case-insensitive. Can be assembly name (e.g. hg38),
            (part of a) scientific name (Danio rerio) or taxonomy id (722).

        Yields
        ------
        tuples with name and metadata
        """
        term = str(term)
        genomes = self.genomes
        if safe(term) in genomes:
            yield self._genome_info_tuple(term)

        elif is_number(term):
            for name in genomes:
                if self._search_taxids(genome=genomes[name], term=term):
                    yield self._genome_info_tuple(name)

        else:
            term = safe(term).lower()
            for name in genomes:
                if self._search_descriptions(genome=genomes[name], term=term):
                    yield self._genome_info_tuple(name)


register_provider = ProviderBase.register_provider


@register_provider("Ensembl")
class EnsemblProvider(ProviderBase):
    """
    Ensembl genome provider.

    Will search both ensembl.org as well as ensemblgenomes.org.
    The bacteria division is not yet supported.
    """

    rest_url = "http://rest.ensembl.org/"

    def __init__(self):
        # Necessary for bucketcache, otherwise methods with identical names
        # from different classes will use the same cache :-O!
        self.name = "Ensembl"
        # Populate on init, so that methods can be cached
        self.genomes = self._get_genomes()
        self.accession_fields = ["assembly_accession"]
        self.taxid_fields = ["taxonomy_id"]
        self.description_fields = [
            "name",
            "scientific_name",
            "url_name",
            "display_name",
        ]
        self.version = None

    @cached(method=True)
    def _request_json(self, ext):
        """Make a REST request and return as json."""
        if self.rest_url.endswith("/") and ext.startswith("/"):
            ext = ext[1:]

        r = requests.get(
            self.rest_url + ext, headers={"Content-Type": "application/json"}
        )

        if not r.ok:
            r.raise_for_status()

        return r.json()

    @cached(method=True)
    def _get_genomes(self):
        sys.stderr.write("Downloading assembly summaries from Ensembl\n")

        genomes = {}
        divisions = self._request_json("info/divisions?")
        for division in divisions:
            if division == "EnsemblBacteria":
                continue
            division_genomes = self._request_json(
                "info/genomes/division/{}?".format(division)
            )
            for genome in division_genomes:
                genomes[safe(genome["assembly_name"])] = genome
        return genomes

    def list_install_options(self, name=None):
        """List Ensembl specific install options"""

        provider_specific_options = {
            "toplevel": {
                "long": "toplevel",
                "help": "always download toplevel-genome",
                "flag_value": True,
            },
            "version": {
                "long": "version",
                "help": "select release version",
                "type": int,
                "default": None,
            },
        }

        return provider_specific_options

    def _genome_info_tuple(self, name):
        """tuple with assembly metadata"""
        genome = self.genomes[name]
        accession = self.assembly_accession(genome)
        taxid = self.genome_taxid(genome)
        taxid = str(taxid) if taxid != 0 else "na"

        return (
            name,
            accession,
            genome.get("scientific_name", "na"),
            taxid,
            genome.get("genebuild", "na"),
        )

    def get_version(self, ftp_site):
        """Retrieve current version from Ensembl FTP.
        """
        with urlopen(ftp_site + "/current_README") as response:
            p = re.compile(r"Ensembl (Genomes|Release) (\d+)")
            m = p.search(response.read().decode())
        if m:
            version = m.group(2)
            sys.stderr.write("Using version {}\n".format(version))
            self.version = version
            return version

    def get_genome_download_link(self, name, mask="soft", **kwargs):
        """
        Return Ensembl http or ftp link to the genome sequence

        Parameters
        ----------
        name : str
            Genome name. Current implementation will fail if exact
            name is not found.

        mask : str , optional
            Masking level. Options: soft, hard or none. Default is soft.

        Returns
        ------
        str with the http/ftp download link.
        """
        genome = self.genomes[safe(name)]

        # parse the division
        division = genome["division"].lower().replace("ensembl", "")
        if division == "bacteria":
            raise NotImplementedError("bacteria from ensembl not yet supported")

        ftp_site = "ftp://ftp.ensemblgenomes.org/pub"
        if division == "vertebrates":
            ftp_site = "http://ftp.ensembl.org/pub"

        # Ensembl release version
        version = self.version
        if kwargs.get("version"):
            version = kwargs.get("version")
        elif not version:
            version = self.get_version(ftp_site)

        # division dependent url format
        ftp_dir = "{}/release-{}/fasta/{}/dna".format(
            division, version, genome["url_name"].lower()
        )
        if division == "vertebrates":
            ftp_dir = "release-{}/fasta/{}/dna".format(
                version, genome["url_name"].lower()
            )
        url = f"{ftp_site}/{ftp_dir}"

        # masking and assembly level
        def get_url(level="toplevel"):
            masks = {"soft": "dna_sm.{}", "hard": "dna_rm.{}", "none": "dna.{}"}
            pattern = masks[mask].format(level)

            asm_url = "{}/{}.{}.{}.fa.gz".format(
                url,
                genome["url_name"].capitalize(),
                re.sub(r"\.p\d+$", "", safe(genome["assembly_name"])),
                pattern,
            )
            return asm_url

        # try to get the (much smaller) primary assembly,
        # unless specified otherwise
        link = get_url("primary_assembly")
        if kwargs.get("toplevel") or not check_url(link):
            link = get_url()

        if check_url(link):
            return link

        raise GenomeDownloadError(
            f"Could not download genome {name} from {self.name}.\n"
            "URL is broken. Select another genome or provider.\n"
            f"Broken URL: {link}"
        )

    def get_annotation_download_link(self, name, **kwargs):
        """
        Parse and test the link to the Ensembl annotation file.

        Parameters
        ----------
        name : str
            Genome name
        kwargs: dict , optional:
            Provider specific options.

            version : int , optional
                Ensembl version. By default the latest version is used.
        """
        genome = self.genomes[safe(name)]
        division = genome["division"].lower().replace("ensembl", "")

        ftp_site = f"ftp://ftp.ensemblgenomes.org/pub"
        if division == "vertebrates":
            ftp_site = "http://ftp.ensembl.org/pub"

        # Ensembl release version
        version = self.version
        if kwargs.get("version"):
            version = kwargs.get("version")
        elif not version:
            version = self.get_version(ftp_site)

        if division != "vertebrates":
            ftp_site += f"/{division}"

        # Get the GTF URL
        base_url = ftp_site + "/release-{}/gtf/{}/{}.{}.{}.gtf.gz"
        safe_name = re.sub(r"\.p\d+$", "", name)
        link = base_url.format(
            version,
            genome["url_name"].lower(),
            genome["url_name"].capitalize(),
            safe_name,
            version,
        )

        if check_url(link):
            return link


@register_provider("UCSC")
class UcscProvider(ProviderBase):
    """
    UCSC genome provider.

    The UCSC API REST server is used to search and list genomes.
    """

    base_url = "http://hgdownload.soe.ucsc.edu/goldenPath"
    ucsc_url = base_url + "/{0}/bigZips/chromFa.tar.gz"
    ucsc_url_masked = base_url + "/{0}/bigZips/chromFaMasked.tar.gz"
    alt_ucsc_url = base_url + "/{0}/bigZips/{0}.fa.gz"
    alt_ucsc_url_masked = base_url + "/{0}/bigZips/{0}.fa.masked.gz"
    rest_url = "http://api.genome.ucsc.edu/list/ucscGenomes"

    def __init__(self):
        # Necessary for bucketcache, otherwise methods with identical names
        # from different classes will use the same cache :-O!
        self.name = "UCSC"
        # Populate on init, so that methods can be cached
        self.genomes = self._get_genomes()
        self.accession_fields = []
        self.taxid_fields = ["taxId"]
        self.description_fields = ["description", "scientificName"]

    @cached(method=True)
    def _get_genomes(self):
        sys.stderr.write("Downloading assembly summaries from UCSC\n")

        r = requests.get(self.rest_url, headers={"Content-Type": "application/json"})
        if not r.ok:
            r.raise_for_status()
        ucsc_json = r.json()
        genomes = ucsc_json["ucscGenomes"]
        return genomes

    # staticmethod, but the decorators cannot be combined
    @cached(method=True)
    def assembly_accession(self, genome):
        """Return the assembly accession (GCA_*) for a genome.

        UCSC does not server the assembly accession through the REST API.
        Therefore, the readme.html is scanned for a GCA assembly id. If it is
        not found, the linked NCBI assembly page will be checked. Especially
        for older genome builds, the GCA will not be present, in which case
        "na" will be returned.

        Parameters
        ----------
        genome : dict
            provider metadata dict of a genome.

        Returns
        ------
        str
            Assembly accession.
        """
        ucsc_url = "https://hgdownload.soe.ucsc.edu/" + genome["htmlPath"]

        p = re.compile(r"GCA_\d+\.\d+")
        p_ncbi = re.compile(r"https?://www.ncbi.nlm.nih.gov/assembly/\d+")
        try:
            text = read_url(ucsc_url)
        except UnicodeDecodeError:
            return "UnicodeDecodeError"
        m = p.search(text)
        # Default, if not found. This matches NCBI, which will also return na.
        gca = "na"
        if m:
            # Get the GCA from the html
            gca = m.group(0)
        else:
            # Search for an assembly link at NCBI
            m = p_ncbi.search(text)
            if m:
                ncbi_url = m.group(0)
                text = read_url(ncbi_url)
                # We need to select the line that contains the assembly accession.
                # The page will potentially contain many more links to newer assemblies
                lines = text.split("\n")
                text = "\n".join(
                    [line for line in lines if "RefSeq assembly accession:" in line]
                )
                m = p.search(text)
                if m:
                    gca = m.group(0)
        return gca

    def _genome_info_tuple(self, name):
        """tuple with assembly metadata"""
        genome = self.genomes[name]
        accession = self.assembly_accession(genome)
        taxid = self.genome_taxid(genome)
        taxid = str(taxid) if taxid != 0 else "na"

        return (
            name,
            accession,
            genome.get("scientificName", "na"),
            taxid,
            genome.get("description", "na"),
        )

    def get_genome_download_link(self, name, mask="soft", **kwargs):
        """
        Return UCSC http link to genome sequence

        Parameters
        ----------
        name : str
            Genome name. Current implementation will fail if exact
            name is not found.

        mask : str , optional
            Masking level. Options: soft, hard or none. Default is soft.

        Returns
        ------
        str with the http/ftp download link.
        """
        # soft masked genomes. can be unmasked in _post _process_download
        urls = [self.ucsc_url, self.alt_ucsc_url]
        if mask == "hard":
            urls = [self.ucsc_url_masked, self.alt_ucsc_url_masked]

        for genome_url in urls:
            link = genome_url.format(name)

            if check_url(link):
                return link

        raise GenomeDownloadError(
            f"Could not download genome {name} from {self.name}.\n"
            "URLs are broken. Select another genome or provider.\n"
            f"Broken URLs: {', '.join([url.format(name) for url in urls])}"
        )

    @staticmethod
    def _post_process_download(name, localname, out_dir, mask="soft"):
        """
        Unmask a softmasked genome if required

        Parameters
        ----------
        name : str
            unused for the UCSC function

        localname : str
            Custom name for your genome

        out_dir : str
            Output directory

        mask : str , optional
            masking level: soft/hard/none, default=soft
        """
        del name
        if mask != "none":
            return

        sys.stderr.write("\nUCSC genomes are softmasked by default. Unmasking...\n")

        fa = os.path.join(out_dir, f"{localname}.fa")
        old_fa = os.path.join(out_dir, f"old_{localname}.fa")
        os.rename(fa, old_fa)
        with open(old_fa) as old, open(fa, "w") as new:
            for line in old:
                if line.startswith(">"):
                    new.write(line)
                else:
                    new.write(line.upper())

    def get_annotation_download_link(self, name, **kwargs):
        """
        Parse and test the link to the UCSC annotation file.

        Will check UCSC, Ensembl and RefSeq annotation, respectively.

        Parameters
        ----------
        name : str
            Genome name
        """
        ucsc_gene_url = f"http://hgdownload.cse.ucsc.edu/goldenPath/{name}/database/"
        annot_files = ["knownGene.txt.gz", "ensGene.txt.gz", "refGene.txt.gz"]

        for file in annot_files:
            link = ucsc_gene_url + file
            if check_url(link):
                return link


@register_provider("NCBI")
class NcbiProvider(ProviderBase):
    """
    NCBI genome provider.

    Uses the assembly reports page to search and list genomes.
    """

    assembly_url = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/"

    def __init__(self):
        # Necessary for bucketcache, otherwise methods with identical names
        # from different classes will use the same cache :-O!
        self.name = "NCBI"
        # Populate on init, so that methods can be cached
        self.genomes = self._get_genomes()
        self.accession_fields = ["assembly_accession", "gbrs_paired_asm"]
        self.taxid_fields = ["species_taxid", "taxid"]
        self.description_fields = ["submitter", "organism_name"]

    @cached(method=True)
    def _get_genomes(self):
        """Parse genomes from assembly summary txt files."""
        sys.stderr.write(
            "Downloading assembly summaries from NCBI, this will take a while...\n"
        )

        genomes = {}
        # order is important as asm_name can repeat (overwriting the older name)
        names = [
            "assembly_summary_refseq_historical.txt",
            "assembly_summary_genbank.txt",
            "assembly_summary_refseq.txt",
        ]
        for fname in names:
            urlcleanup()
            with urlopen(os.path.join(self.assembly_url, fname)) as response:
                lines = response.read().decode("utf-8").splitlines()
            header = lines[1].strip("# ").split("\t")
            for line in lines[2:]:
                vals = line.strip("# ").split("\t")
                # overwrites older asm_names
                genomes[safe(vals[15])] = dict(zip(header, vals))
        return genomes

    def _genome_info_tuple(self, name):
        """tuple with assembly metadata"""
        genome = self.genomes[name]
        accession = self.assembly_accession(genome)
        taxid = self.genome_taxid(genome)
        taxid = str(taxid) if taxid != 0 else "na"

        return (
            name,
            accession,
            genome.get("organism_name", "na"),
            taxid,
            genome.get("submitter", "na"),
        )

    def get_genome_download_link(self, name, mask="soft", **kwargs):
        """
        Return NCBI ftp link to top-level genome sequence

        Parameters
        ----------
        name : str
            Genome name. Current implementation will fail if exact
            name is not found.

        mask : str , optional
            Masking level. Options: soft, hard or none. Default is soft.

        Returns
        ------
        str with the http/ftp download link.
        """
        genome = self.genomes[safe(name)]

        # only soft masked genomes available. can be (un)masked in _post _process_download
        link = genome["ftp_path"]
        link = link.replace("ftp://", "https://")
        link += "/" + link.split("/")[-1] + "_genomic.fna.gz"

        if check_url(link):
            return link

        raise GenomeDownloadError(
            f"Could not download genome {name} from {self.name}.\n"
            "URL is broken. Select another genome or provider.\n"
            f"Broken URL: {link}"
        )

    def _post_process_download(self, name, localname, out_dir, mask="soft"):
        """
        Replace accessions with sequence names in fasta file.

        Applies masking.

        Parameters
        ----------
        name : str
            NCBI genome name

        localname : str
            Custom name for your genome

        out_dir : str
            Output directory

        mask : str , optional
            masking level: soft/hard/none, default=soft
        """
        # Create mapping of accessions to names
        genome = self.genomes[safe(name)]
        url = genome["ftp_path"]
        url += f"/{url.split('/')[-1]}_assembly_report.txt"
        url = url.replace("ftp://", "https://")

        tr = {}
        urlcleanup()
        with urlopen(url) as response:
            for line in response.read().decode("utf-8").splitlines():
                if line.startswith("#"):
                    continue
                vals = line.strip().split("\t")
                tr[vals[6]] = vals[0]

        # mask sequence if required
        if mask == "soft":

            def mask_cmd(txt):
                return txt

        elif mask == "hard":
            sys.stderr.write(
                "\nNCBI genomes are softmasked by default. Hard masking...\n"
            )

            def mask_cmd(txt):
                return re.sub("[actg]", "N", txt)

        else:
            sys.stderr.write("\nNCBI genomes are softmasked by default. Unmasking...\n")

            def mask_cmd(txt):
                return txt.upper()

        # apply mapping and masking
        fa = os.path.join(out_dir, f"{localname}.fa")
        old_fa = os.path.join(out_dir, f"old_{localname}.fa")
        os.rename(fa, old_fa)
        with open(old_fa) as old, open(fa, "w") as new:
            for line in old:
                if line.startswith(">"):
                    desc = line.strip()[1:]
                    name = desc.split(" ")[0]
                    new.write(">{} {}\n".format(tr.get(name, name), desc))
                else:
                    new.write(mask_cmd(line))

    def get_annotation_download_link(self, name, **kwargs):
        """
        Parse and test the link to the NCBI annotation file.

        Parameters
        ----------
        name : str
            Genome name
        """
        genome = self.genomes[safe(name)]
        link = genome["ftp_path"]
        link = link.replace("ftp://", "https://")
        link += "/" + link.split("/")[-1] + "_genomic.gff.gz"

        if check_url(link):
            return link


@register_provider("URL")
class UrlProvider(ProviderBase):
    """
    URL genome provider.

    Simply download a genome directly through an url.
    """

    def __init__(self):
        self.name = "URL"
        self.genomes = {}

    def genome_taxid(self, genome):
        return "na"

    def assembly_accession(self, genome):
        return "na"

    def search(self, term):
        """return an empty generator,
        same as if no genomes were found at the other providers"""
        yield from ()

    def list_install_options(self, name=None):
        """List URL specific install options"""

        provider_specific_options = {
            "to_annotation": {
                "long": "to-annotation",
                "help": "link to the annotation file, required if this is not in the same directory as the fasta file",
                "default": None,
            },
        }

        return provider_specific_options

    def get_genome_download_link(self, url, mask=None, **kwargs):
        return url

    def get_annotation_download_link(self, name, **kwargs):
        """
        check if the linked annotation file is of a supported file type (gtf/gff3/bed)
        """
        link = kwargs.get("to_annotation")
        if link:
            ext = get_file_info(link)[0]
            if ext not in [".gtf", ".gff", ".gff3", ".bed"]:
                raise TypeError(
                    "Only (gzipped) gtf, gff and bed files are supported.\n"
                )

            if check_url(link):
                return link

    @staticmethod
    def search_url_for_annotation(url):
        """Attempts to find a gtf or gff3 file in the same location as the genome url"""
        urldir = os.path.dirname(url)
        sys.stderr.write(
            "You have requested gene annotation to be downloaded.\n"
            "Genomepy will check the remote directory:\n"
            f"{urldir} for annotation files...\n"
        )

        # try to find a GTF or GFF3 file
        name = get_localname(url)
        with urlopen(urldir) as f:
            for urlline in f.readlines():
                urlstr = str(urlline)
                if any(
                    substring in urlstr.lower() for substring in [".gtf", name + ".gff"]
                ):
                    break

        # retrieve the filename from the HTML line
        fname = ""
        for split in re.split('>|<|><|/|"', urlstr):
            if split.lower().endswith(
                (
                    ".gtf",
                    ".gtf.gz",
                    name + ".gff",
                    name + ".gff.gz",
                    name + ".gff3",
                    name + ".gff3.gz",
                )
            ):
                fname = split
                break
        else:
            raise FileNotFoundError(
                "Could not parse the remote directory. "
                "Please supply a URL using --url-to-annotation.\n"
            )

        # set variables for downloading
        link = urldir + "/" + fname

        if check_url(link):
            return link

    def download_annotation(self, url, genomes_dir, localname=None, **kwargs):
        """
        Attempts to download a gtf or gff3 file from the same location as the genome url

        Parameters
        ----------
        url : str
            url of where to download genome from

        genomes_dir : str
            Directory to install annotation

        localname : str , optional
            Custom name for your genome

        kwargs: dict , optional:
            Provider specific options.

            to_annotation : str , optional
                url to annotation file (only required if this not located in the same directory as the fasta)
        """
        name = get_localname(url)
        localname = get_localname(name, localname)

        if kwargs.get("to_annotation"):
            link = self.get_annotation_download_link(None, **kwargs)
        else:
            link = self.search_url_for_annotation(url)

        self.attempt_and_report(name, localname, link, genomes_dir)
