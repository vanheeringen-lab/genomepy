import os
from time import sleep

from loguru import logger

from genomepy.caching import cache
from genomepy.exceptions import GenomeDownloadError
from genomepy.files import update_readme
from genomepy.online import check_url, connect_ftp_link
from genomepy.providers.base import BaseProvider, download_annotation
from genomepy.providers.ncbi import download_assembly_report
from genomepy.providers.ucsc import UcscProvider
from genomepy.utils import get_genomes_dir, get_localname


class GencodeProvider(BaseProvider):
    """
    GENCODE genome provider.

    GENCODE sports superb annotations for human and mouse with UCSC-style chromosome names.
    Genomes on this provider are unmasked, so we use the UCSC genomes instead.

    Note: this combination lacks scaffolds and alternate haplotype sequences
    (their names don't match between GENCODE and UCSC).
    """

    name = "GENCODE"
    accession_fields = ["assembly_accession"]
    taxid_fields = ["taxonomy_id"]
    description_fields = [
        "species",
        "other_info",
        "text_search",
    ]
    _cli_install_options = {}
    _ftp_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode"

    def __init__(self):
        self._provider_status()
        # Populate on init, so that methods can be cached
        self.genomes = _get_genomes(self._ftp_link)
        self.ucsc = UcscProvider()
        self.gencode2ucsc = get_gencode2ucsc(self.genomes)
        self._update_genomes()

    @staticmethod
    def ping():
        """Can the provider be reached?"""
        return bool(check_url("ftp.ebi.ac.uk/pub/databases/gencode"))

    def _genome_info_tuple(self, name):
        """tuple with assembly metadata"""
        accession = self.genomes[name]["assembly_accession"]
        taxid = self.genomes[name]["taxonomy_id"]
        annotations = True
        species = self.genomes[name]["species"]
        other = self.genomes[name]["other_info"]
        return name, accession, taxid, annotations, species, other

    def _update_genomes(self):
        """add assembly accession and other information to the genomes dict"""
        for name in self.genomes:
            self.genomes[name]["other_info"] = "GENCODE annotation + UCSC genome"
            # makes the genome findable with UCSC names
            ucsc_name = self.gencode2ucsc[name]
            self.genomes[name]["text_search"] += f" {ucsc_name}"
            # add the UCSC accession ID
            ucsc_acc = self.ucsc.genomes[ucsc_name]["assembly_accession"]
            self.genomes[name]["assembly_accession"] = ucsc_acc

    def get_genome_download_link(self, name: str, mask: str = "soft", **kwargs) -> str:
        """
        Return UCSC http link to genome sequence

        Parameters
        ----------
        name : str
            Genome name.

        mask : str , optional
            Masking level. Options: soft, hard or none. Default is soft.

        Returns
        ------
        str
            http/ftp link tp genome.
        """
        ucsc_name = self.gencode2ucsc[name]
        return self.ucsc.get_genome_download_link(ucsc_name, mask, **kwargs)

    def download_genome(
        self,
        name: str,
        genomes_dir: str = None,
        localname: str = None,
        mask: str = "soft",
        **kwargs,
    ):
        """
        Download genomes from UCSC, as the GENCODE genomes aren't masked.
        Contigs between the UCSC genome and GENCODE annotations match.
        """
        ucsc_name = self.gencode2ucsc[name]
        self.ucsc.name = "GENCODE"  # for logging & readme
        self.ucsc.download_genome(ucsc_name, genomes_dir, localname, mask, **kwargs)
        self.ucsc.name = "UCSC"

    def get_annotation_download_links(self, name, **kwargs):
        """
        Retrieve functioning gene annotation download link(s).

        Parameters
        ----------
        name : str
            genome name

        Returns
        -------
        list
            http/ftp link(s)
        """
        return self.genomes[name]["annotations"]

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
        name = self._check_name(name)
        link = self.get_annotation_download_link(name, **kwargs)

        localname = get_localname(name, localname)
        genomes_dir = get_genomes_dir(genomes_dir, check_exist=False)

        logger.info(f"Downloading annotation from {self.name}. Target URL: {link}...")
        try:
            # download exact assembly report to rename the scaffolds
            acc = self.assembly_accession(name)
            fname = os.path.join(genomes_dir, localname, "assembly_report.txt")
            download_assembly_report(acc, fname)

            download_annotation(genomes_dir, link, localname)
            logger.info("Annotation download successful")
        except Exception as e:
            raise GenomeDownloadError(
                f"An error occured while installing the gene annotation for {name} from {self.name}.\n"
                "If you think the annotation should be there, please file a bug report at: "
                "https://github.com/vanheeringen-lab/genomepy/issues\n\n"
                f"Error: {e.args[0]}"
            )

        # Add annotation URL to readme
        readme = os.path.join(genomes_dir, localname, "README.txt")
        update_readme(readme, updated_metadata={"annotation url": link})


def get_gencode2ucsc(genomes):
    """
    generate conversion dict for gencode's
    Ensembl-style names to UCSC-style names.

    Parameters
    ----------
    genomes : dict
        the provider's genomes dict

    Returns
    -------
    dict
        keys: Ensembl names, values: UCSC names
    """
    # start with historic assemblies
    gencode2ucsc = {
        "GRCh37": "hg19",
        "GRCm38": "mm10",
    }
    # add latest assemblies
    for assembly in genomes:
        if assembly not in gencode2ucsc:
            specie = genomes[assembly]["species"]
            specie = "hg" if specie.startswith("H") else "mm"
            number = "".join([s for s in assembly if s.isdigit()])
            ucsc_assembly = f"{specie}{number}"
            gencode2ucsc[assembly] = ucsc_assembly
    return gencode2ucsc


def get_releases(listing: list, specie: str) -> list:
    """
    return the list of gencode releases found on the site.
    Highest fist.

    Parameters
    ----------
    listing : list
        list of gencode FTP paths, some ending with "release_nn" or "release_Mnn"

    specie : str
        mouse of human. mouse release numbers are prefixed with "M"

    Returns
    ------
    list
        releases sorted from highest to lowest
    """
    # ignore releases without a consistent system.
    releases = [int(ftpdir[-2:]) for ftpdir in listing if ftpdir[-2:].isdigit()]
    releases = sorted([str(n) for n in releases if n > 21], reverse=True)
    if specie == "mouse":
        releases = ["M" + r for r in releases]  # dear gencode: why?
    return releases


def add_grch37(genomes, ftp_link):
    """gencode's GRCh37 is filed with a unique system"""
    # no mouse in genomes yet
    latest_human = f"GRCh{max([int(g[-2:]) for g in genomes])}"
    latest_annot = genomes[latest_human]["annotations"][0]
    release = [r for r in latest_annot.split("/") if "release_" in r][0][-2:]
    genomes["GRCh37"] = {
        "annotations": [
            f"{ftp_link}/Gencode_human/release_{release}/GRCh37_mapping/"
            f"gencode.v{release}lift37.annotation.gtf.gz"
        ],
        "taxonomy_id": 9606,
        "species": "Homo sapiens",
        "text_search": "human",
    }
    return genomes


def _get_genomes(ftp_link):
    """Wrapper to catch FTP issues."""
    try:
        genomes = get_genomes(ftp_link)
    except (TimeoutError, ConnectionRefusedError):
        logger.warning(
            "GENCODE is online, but cannot be reached. "
            "Is FTP working on this device?"
        )
        genomes = {}
    return genomes


@cache
def get_genomes(ftp_link):
    """genomes dict of the latest gencode release of each major assembly."""
    logger.info("Downloading assembly summaries from GENCODE")

    genomes = {}
    species = {"human": "Homo sapiens", "mouse": "Mus musculus"}
    taxid = {"human": 9606, "mouse": 10090}
    sleep(1)  # we just closed a connection with ping()
    ftp, ftp_path = connect_ftp_link(ftp_link, timeout=7)
    for specie in ["human", "mouse"]:
        listing = ftp.nlst(f"{ftp_path}/Gencode_{specie}")
        releases = get_releases(listing, specie)
        for release in releases:
            listing = ftp.nlst(f"{ftp_path}/Gencode_{specie}/release_{release}")
            files = [f.split("/")[-1] for f in listing]
            # drop patch level (UCSC genomes don't have patches either)
            assembly = [f for f in files if "primary_assembly" in f][0].split(".")[0]
            if assembly not in genomes:
                genomes[assembly] = {
                    "annotations": [
                        f"{ftp_link}/Gencode_{specie}/release_{release}/"
                        f"gencode.v{release}.primary_assembly.annotation.gtf.gz"
                    ],
                    "taxonomy_id": taxid[specie],
                    "species": species[specie],
                    "text_search": specie,
                }

        if specie == "human":
            genomes = add_grch37(genomes, ftp_link)

    ftp.quit()
    return genomes
