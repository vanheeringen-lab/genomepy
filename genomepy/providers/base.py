import os
import time
import shutil
import subprocess as sp

from appdirs import user_cache_dir
from bucketcache import Bucket
from loguru import logger
from tempfile import mkdtemp
from typing import Iterator, Union

from genomepy.__about__ import __version__
from genomepy.exceptions import GenomeDownloadError
from genomepy.files import (
    write_readme,
    update_readme,
    gunzip_and_name,
    tar_to_bigfile,
    get_file_info,
)
from genomepy.online import download_file, check_url
from genomepy.utils import (
    get_localname,
    safe,
    lower,
    mkdir_p,
    get_genomes_dir,
    rm_rf,
)


# Store the output of slow commands (marked with @cache and @goldfish_cache) for fast reuse.
# Bucketcache creates a new pickle for each function + set of unique variables,
# to avoid storing self.genomes, ignore=["self"] or use @staticmethod.
my_cache_dir = os.path.join(user_cache_dir("genomepy"), __version__)
if not os.path.exists(my_cache_dir):
    os.makedirs(my_cache_dir)
cache = Bucket(my_cache_dir, days=7)
goldfish_cache = Bucket(my_cache_dir, minutes=10)


class BaseProvider:
    """
    Provider base class.
    """

    # class variables set by child classes:
    name = None
    genomes = {}
    accession_fields = []
    taxid_fields = []
    description_fields = []
    provider_specific_install_options = {}

    def __hash__(self):
        return hash(str(self.__class__))

    @goldfish_cache(ignore=["self", "max_tries"])
    def provider_status(self, url, max_tries=1):
        """check if provider is online (stores results for 10 minutes)"""
        if not check_url(url, max_tries):
            raise ConnectionError(f"{self.name} appears to be offline.\n")

    def check_name(self, name):
        """check if genome name can be found for provider"""
        if not safe(name) in self.genomes:
            raise GenomeDownloadError(
                f"Could not download genome {name} from {self.name}.\n\n"
                "Check for typos or try\n"
                f"  genomepy search {name} -p {self.name}\n"
            )

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
        for name in self.genomes.keys():
            yield self._genome_info_tuple(name)

    def genome_taxid(self, name: str) -> int:
        """
        Return the genome taxonomy ID for a genome.

        Parameters
        ----------
        name: str
            genome name

        Returns
        ------
        int
            Genome Taxonomy identifier
        """
        for field in self.taxid_fields:
            tid = str(self.genomes[name].get(field))
            if tid.isdigit():
                return int(tid)

    def assembly_accession(self, name: str) -> str:
        """
        Return the assembly accession number (GCA* or GCF*) for a genome.

        Parameters
        ----------
        name: str
            genome name

        Returns
        ------
        str
            Assembly accession number
        """
        for field in self.accession_fields:
            accession = str(self.genomes[name].get(field))
            if accession.startswith(("GCA", "GCF")):
                return accession

    def annotation_links(self, name: str) -> list:
        """
        Return available gene annotation links (http/ftp) for a genome

        Parameters
        ----------
        name: str
            genome name

        Returns
        ------
        list
            Gene annotation links
        """
        if "annotations" not in self.genomes[safe(name)]:
            links = self.get_annotation_download_links(name)
            self.genomes[safe(name)]["annotations"] = links
        return self.genomes[safe(name)]["annotations"]

    def get_genome_download_link(self, name, mask="soft", **kwargs):
        raise NotImplementedError()

    def download_genome(
        self,
        name: str,
        genomes_dir: str = None,
        localname: str = None,
        mask: str = "soft",
        **kwargs,
    ):
        """
        Download a (gzipped) genome file to a specific directory

        Parameters
        ----------
        name : str
            Genome / species name

        genomes_dir : str , optional
            Directory to install genome

        localname : str , optional
            Custom name for your genome

        mask: str , optional
            Masking, soft, hard or none (all other strings)
        """
        name = safe(name)
        self.check_name(name)

        link = self.get_genome_download_link(name, mask=mask, **kwargs)

        localname = get_localname(name, localname)
        genomes_dir = get_genomes_dir(genomes_dir, check_exist=False)
        out_dir = os.path.join(genomes_dir, localname)
        mkdir_p(out_dir)

        logger.info(f"Downloading genome from {self.name}. Target URL: {link}...")

        # download to tmp dir. Move genome on completion.
        # tmp dir is in genome_dir to prevent moving the genome between disks
        tmp_dir = mkdtemp(dir=out_dir)
        fname = os.path.join(tmp_dir, f"{localname}.fa")

        download_file(link, fname)
        logger.info("Genome download successful, starting post processing...")

        # unzip genome
        if link.endswith(".tar.gz"):
            tar_to_bigfile(fname, fname)
        elif link.endswith(".gz"):
            os.rename(fname, fname + ".gz")
            gunzip_and_name(fname + ".gz")

        # process genome (e.g. masking)
        if hasattr(self, "_post_process_download"):
            self._post_process_download(
                name=name, localname=localname, out_dir=tmp_dir, mask=mask
            )

        # transfer the genome from the tmpdir to the genome_dir
        src = fname
        dst = os.path.join(genomes_dir, localname, os.path.basename(fname))
        shutil.move(src, dst)
        rm_rf(tmp_dir)

        asm_report = os.path.join(out_dir, "assembly_report.txt")
        asm_acc = self.assembly_accession(name)
        if hasattr(self, "download_assembly_report") and asm_acc:
            self.download_assembly_report(asm_acc, asm_report)

        logger.info("name: {}".format(name))
        logger.info("local name: {}".format(localname))
        logger.info("fasta: {}".format(dst))

        # Create readme with information
        readme = os.path.join(genomes_dir, localname, "README.txt")
        metadata = {
            "name": localname,
            "provider": self.name,
            "original name": name,
            "original filename": os.path.split(link)[-1],
            "assembly_accession": asm_acc if asm_acc else "na",
            "tax_id": self.genome_taxid(name) if self.genome_taxid(name) else "na",
            "mask": mask,
            "genome url": link,
            "annotation url": "na",
            "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        write_readme(readme, metadata)

    def get_annotation_download_links(self, name, **kwargs):
        raise NotImplementedError()

    def get_annotation_download_link(self, name: str, **kwargs) -> str:
        """select a functional annotation download link from a list of links"""
        links = self.annotation_links(name)
        if links:
            return links[0]
        logger.error(f"No gene annotations found for {name} on {self.name}")
        raise FileNotFoundError

    @staticmethod
    def _download_annotation(genomes_dir, annot_url, localname):
        """download annotation file, convert to intermediate file and generate output files"""

        # create output directory if missing
        out_dir = os.path.join(genomes_dir, localname)
        mkdir_p(out_dir)

        # download to tmp dir. Move genome on completion.
        # tmp dir is in genome_dir to prevent moving the genome between disks
        tmp_dir = mkdtemp(dir=out_dir)
        ext, gz = get_file_info(annot_url)
        annot_file = os.path.join(tmp_dir, localname + ".annotation" + ext)
        download_file(annot_url, annot_file)

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
            cmd = "gtfToGenePred -ignoreGroupsWithoutExons {0} {1}"
        elif "txt" in ext:
            # UCSC annotations only
            with open(annot_file) as f:
                cols = f.readline().split("\t")

            # extract the genePred format columns
            start_col = 1
            for i, col in enumerate(cols):
                if col in ["+", "-"]:
                    start_col = i - 1
                    break
            end_col = start_col + 10
            cmd = (
                f"""cat {{0}} | cut -f {start_col}-{end_col} | """
                # knownGene.txt.gz has spotty fields, this replaces non-integer fields with zeroes
                + """awk 'BEGIN {{FS=OFS="\t"}} !($11 ~ /^[0-9]+$/) {{$11="0"}}1' > {1}"""
            )
        else:
            raise TypeError(f"file type extension {ext} not recognized!")

        sp.check_call(cmd.format(annot_file, pred_file), shell=True)

        # generate gzipped gtf file (if required)
        gtf_file = annot_file.replace(ext, ".gtf")
        if "gtf" not in ext:
            cmd = "genePredToGtf -source=genomepy file {0} {1}"
            sp.check_call(cmd.format(pred_file, gtf_file), shell=True)

        # generate gzipped bed file (if required)
        bed_file = annot_file.replace(ext, ".bed")
        if "bed" not in ext:
            cmd = "genePredToBed {0} {1}"
            sp.check_call(cmd.format(pred_file, bed_file), shell=True)

        # transfer the files from the tmpdir to the genome_dir
        for f in [gtf_file, bed_file]:
            src = f
            dst = os.path.join(out_dir, os.path.basename(f))
            shutil.move(src, dst)
        rm_rf(tmp_dir)

    def download_annotation(self, name, genomes_dir=None, localname=None, **kwargs):
        """
        Download annotation file to to a specific directory

        Parameters
        ----------
        name : str
            Genome / species name

        genomes_dir : str , optional
            Directory to install annotation

        localname : str , optional
            Custom name for your genome
        """
        self.check_name(name)

        link = self.get_annotation_download_link(name, **kwargs)

        localname = get_localname(name, localname)
        genomes_dir = get_genomes_dir(genomes_dir, check_exist=False)

        logger.info(f"Downloading annotation from {self.name}. Target URL: {link}...")
        try:
            self._download_annotation(genomes_dir, link, localname)
            logger.info("Annotation download successful")
        except Exception:
            raise GenomeDownloadError(
                f"\nCould not download annotation for {name} from {self.name}\n"
                "If you think the annotation should be there, please file a bug report at:\n"
                "https://github.com/vanheeringen-lab/genomepy/issues\n"
            )

        # Add annotation URL to readme
        readme = os.path.join(genomes_dir, localname, "README.txt")
        update_readme(readme, updated_metadata={"annotation url": link})

    def _search_text(self, term: str) -> Iterator[str]:
        """check if search term is found in the provider's genome name or description field(s)"""
        for name, metadata in self.genomes.items():
            if term in lower(name) or any(
                [term in lower(metadata[f]) for f in self.description_fields]
            ):
                yield name

    def _search_accession(self, term: str) -> Iterator[str]:
        """check if search term is found in the provider's accession field(s)"""
        # cut off prefix (GCA_/GCF_) and suffix (version numbers, e.g. '.3')
        term = term[4:].split(".")[0]
        for name, metadata in self.genomes.items():
            if any([term in str(metadata[f]) for f in self.accession_fields]):
                yield name

    def _search_taxonomy(self, term: str) -> Iterator[str]:
        """check if search term matches to any of the provider's taxonomy field(s)"""
        for name, metadata in self.genomes.items():
            if any([term == lower(metadata[f]) for f in self.taxid_fields]):
                yield name

    def search(self, term: Union[str, int]):
        """
        Search for term in genome names, descriptions and taxonomy ID.

        The search is case-insensitive.

        Parameters
        ----------
        term : str, int
            Search term, case-insensitive.
            Can be (part of) an assembly name (e.g. hg38),
            scientific name (Danio rerio) or assembly
            accession (GCA_000146045/GCF_...),
            or an exact taxonomy id (7227).

        Yields
        ------
        tuples with name and metadata
        """
        term = lower(term)

        search_function = self._search_text
        if term.startswith(("gca_", "gcf_")):
            search_function = self._search_accession
        if term.isdigit():
            search_function = self._search_taxonomy

        for name in search_function(term):
            yield self._genome_info_tuple(name)
