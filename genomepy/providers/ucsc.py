import os
import re
import subprocess as sp
import urllib.error
from tempfile import mkdtemp
from typing import Generator, Iterator, List

import mysql.connector
import pandas as pd
import requests
from loguru import logger

from genomepy.caching import cache
from genomepy.exceptions import GenomeDownloadError
from genomepy.files import update_readme
from genomepy.online import check_url, read_url
from genomepy.providers.base import BaseProvider
from genomepy.providers.ncbi import NcbiProvider
from genomepy.utils import get_genomes_dir, get_localname, lower, mkdir_p, rm_rf

# order determines which annotation genomepy will attempt to install
# for more info, see http://genome.ucsc.edu/FAQ/FAQgenes.html
ANNOTATIONS = ["ncbiRefSeq", "refGene", "ensGene", "knownGene"]


class UcscProvider(BaseProvider):
    """
    UCSC genome provider.

    The UCSC API REST server is used to search and list genomes.
    The UCSC MySQL database is used to find metadata and annotations.
    """

    name = "UCSC"
    accession_fields = ["assembly_accession"]
    taxid_fields = ["taxId"]
    description_fields = ["description", "scientificName"]
    _cli_install_options = {
        "ucsc_annotation_type": {
            "long": "annotation",
            "help": "specify annotation to download: "
            f"{', '.join(ANNOTATIONS)} (case-insensitive)",
            "default": None,
        },
    }
    _url = "http://hgdownload.soe.ucsc.edu/goldenPath"

    def __init__(self):
        self._provider_status()
        # Populate on init, so that methods can be cached
        self.genomes = get_genomes("http://api.genome.ucsc.edu/list/ucscGenomes")

    @staticmethod
    def ping():
        """Can the provider be reached?"""
        return bool(check_url("http://hgdownload.soe.ucsc.edu/goldenPath"))

    def _search_accession(self, term: str) -> Iterator[str]:
        """
        UCSC does not always store assembly accessions.
        If no hits were found, search accession on NCBI (most genomes + stable accession IDs),
        then uses those search results for a UCSC text search.

        Parameters
        ----------
        term : str
            Assembly accession, `GCA_`/`GCF_`

        Yields
        ------
        genome names
        """
        # cut off prefix (GCA_/GCF_) and suffix (version numbers, e.g. '.3')
        term = term[4:].split(".")[0]
        hits = 0
        for name, metadata in self.genomes.items():
            if any([term in str(metadata[f]) for f in self.accession_fields]):
                hits += 1
                yield name

        # search NCBI only if we found no local hits
        if hits == 0:
            return self._search_accession_ncbi(term)

    def _search_accession_ncbi(self, term: str) -> Iterator[str]:
        """
        search NCBI (most genomes + stable accession IDs),
        then uses the NCBI accession search results for a UCSC text search.

        Parameters
        ----------
        term : str
            Assembly accession, `GCA_`/`GCF_`

        Yields
        ------
        genome names
        """
        # NCBI provides a consistent assembly accession. This can be used to
        # retrieve the species, and then search for that.
        p = NcbiProvider()
        ncbi_genomes = list(p._search_accession(term))

        # remove superstrings (keep GRCh38, not GRCh38.p1 to GRCh38.p13)
        unique_ncbi_genomes = []
        for i in ncbi_genomes:
            if sum([j in i for j in ncbi_genomes]) == 1:
                unique_ncbi_genomes.append(i)

        # add NCBI organism names to search terms
        organism_names = [
            p.genomes[name]["organism_name"] for name in unique_ncbi_genomes
        ]
        terms = list(set(unique_ncbi_genomes + organism_names))

        # search with NCBI results in the given provider
        for name, metadata in self.genomes.items():
            for term in terms:
                term = lower(term)
                if term in lower(name) or any(
                    [term in lower(metadata[f]) for f in self.description_fields]
                ):
                    yield name
                    break  # max one hit per genome

    def assembly_accession(self, name: str) -> str:
        """
        Return the assembly accession (`GCA_`/`GCF_`) for a genome.

        Some accession IDs can be retrieved from the UCSC MySQL hgFixed database.
        For others, the accession IDs can sometimes be scraped from the readme.html.
        If not, any linked NCBI assembly pages can also be scraped.

        Parameters
        ----------
        name: str
            genome name

        Returns
        ------
        str
            Assembly accession.
        """
        acc = self.genomes[name]["assembly_accession"]
        if acc:
            return acc

        acc = scrape_accession(self.genomes[name]["htmlPath"])
        if acc:
            return acc

        return "na"

    def annotation_links(self, name, **kwargs) -> List[str]:
        """
        Return a sorted list of available gene annotation types for a genome

        Parameters
        ----------
        name: str
            genome name

        Returns
        ------
        list
            Gene annotation types
        """
        available = self.get_annotation_download_links(name)
        annotations_found = [a for a in ANNOTATIONS if a in available]
        return annotations_found

    def _genome_info_tuple(self, name):
        """tuple with assembly metadata"""
        accession = self.assembly_accession(name)
        taxid = self.genome_taxid(name)
        annotations = [a in self.annotation_links(name) for a in ANNOTATIONS]
        species = self.genomes[name].get("scientificName")
        other = self.genomes[name].get("description")
        return name, accession, taxid, annotations, species, other

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
        ucsc_url = self._url + "/{0}/bigZips/chromFa.tar.gz"
        ucsc_url_masked = self._url + "/{0}/bigZips/chromFaMasked.tar.gz"
        alt_ucsc_url = self._url + "/{0}/bigZips/{0}.fa.gz"
        alt_ucsc_url_masked = self._url + "/{0}/bigZips/{0}.fa.masked.gz"

        # soft masked genomes. can be unmasked in _post _process_download
        urls = [ucsc_url, alt_ucsc_url]
        if mask == "hard":
            urls = [ucsc_url_masked, alt_ucsc_url_masked]

        for genome_url in urls:
            link = genome_url.format(name)

            if check_url(link, 2):
                return link

        raise GenomeDownloadError(
            f"Could not download genome {name} from {self.name}.\n"
            "URLs are broken. Select another genome or provider.\n"
            f"Broken URLs: {', '.join([url.format(name) for url in urls])}"
        )

    @staticmethod
    def _post_process_download(name, fname, out_dir, mask="soft"):  # noqa
        """
        Unmask a softmasked genome if required

        Parameters
        ----------
        name : str
            unused for the UCSC function

        fname : str
            file path to the genome fasta

        out_dir : str
            unused for the UCSC function

        mask : str , optional
            masking level: soft/hard/none, default=soft
        """
        if mask != "none":
            return

        logger.info("UCSC genomes are softmasked by default. Unmasking...")

        old_fname = os.path.join(
            os.path.dirname(fname), f"original_{os.path.basename(fname)}"
        )
        os.rename(fname, old_fname)
        with open(old_fname) as old, open(fname, "w") as new:
            for line in old:
                if line[0] == ">":
                    new.write(line)
                else:
                    new.write(line.upper())

    def get_annotation_download_links(self, name, **kwargs):
        """
        Return available gene annotation table(s) from the UCSC MySQL database.

        Available tables were retrieved on init.

        Parameters
        ----------
        name : str
            genome name

        Returns
        -------
        list
            annotation types
        """
        return self.genomes[name]["annotations"]

    def get_annotation_download_link(self, name: str, **kwargs) -> str:
        """
        Return an available annotation type.

        Parameters
        ----------
        name : str
            genome name
        **kwargs: dict, optional:
            ucsc_annotation_type : specific annotation type to download.

        Returns
        -------
        str
            http/ftp link

        Raises
        ------
        GenomeDownloadError
            if no functional link was found
        FileNotFoundError
            if the specified annotation type is unavailable
        """
        available = self.annotation_links(name)
        if not available:
            raise GenomeDownloadError(
                f"No gene annotations found for {name} on {self.name}.\n"
                "Check for typos or try\n"
                f"  genomepy search {name} -p {self.name}"
            )
        annot = available

        usr_annot = kwargs.get("ucsc_annotation_type")
        if usr_annot:
            # not all types are available for each genome
            annot = [a for a in available if a.lower() == usr_annot.lower()]
            if not annot:
                raise FileNotFoundError(
                    f"{usr_annot} is not available for {name}. "
                    f"Options: {', '.join(available)}.\n"
                )

        return annot[0]

    def download_annotation(self, name, genomes_dir=None, localname=None, **kwargs):
        """
        Download the UCSC genePred via their MySQL database, and convert to annotations.
        """
        name = self._check_name(name)
        annot = self.get_annotation_download_link(name, **kwargs)

        localname = get_localname(name, localname)
        genomes_dir = get_genomes_dir(genomes_dir, check_exist=False)

        logger.info(f"Downloading the {annot} annotation from the UCSC MySQL database.")
        try:
            download_annotation(name, annot, genomes_dir, localname)
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
        update_readme(
            readme,
            updated_metadata={
                "annotation url": f"UCSC MySQL database: {name}, table: {annot}"
            },
        )

    def head_annotation(self, name, genomes_dir=None, n: int = 5, **kwargs):
        """
        Download the first n genes of each UCSC annotation type.

        The first line of the GTF is printed for review
        (of the gene_name field, for instance).

        Parameters
        ----------
        name : str
            genome name
        genomes_dir : str, optional
            genomes directory to install the annotation in.
        n : int, optional
            download the annotation for n genes.
        kwargs : dict , optional
            annotations : list
                specify which UCSC annotation types to download.
                Downloads all available if left blank.
        """
        name = self._check_name(name)
        genomes_dir = get_genomes_dir(genomes_dir, check_exist=False)

        annotations = kwargs.get("annotations")
        if annotations is None:
            annotations = self.annotation_links(name)

        for annot in annotations:
            if annot not in ANNOTATIONS:
                raise ValueError(f"{annot} not in {ANNOTATIONS}")

            localname = f"{name}_head_{annot}"
            fpath = os.path.join(genomes_dir, localname, f"{localname}.annotation.gtf")
            download_annotation(name, annot, genomes_dir, localname, n=n)

            logger.info(f"{self.name} {annot}")
            with open(fpath) as f:
                for m, line in enumerate(f):
                    if line:
                        print(line.strip())
                        if m + 1 == n:
                            break


@cache
def get_genomes(rest_url):
    logger.info("Downloading assembly summaries from UCSC")

    r = requests.get(rest_url, headers={"Content-Type": "application/json"})
    if not r.ok:
        r.raise_for_status()
    ucsc_json = r.json()
    genomes = ucsc_json["ucscGenomes"]

    for genome in genomes:
        genomes[genome]["assembly_accession"] = None
        genomes[genome]["annotations"] = []
    # add accession IDs (self.assembly_accession will try to fill in the blanks)
    genomes = add_accessions1(genomes)
    genomes = add_accessions2(genomes)
    genomes = add_annotation_links(genomes)

    return genomes


def add_accessions1(genomes: dict) -> dict:
    """
    Use the UCSC MySQL database to obtain accession IDs
    (of the most similar assemblies on other providers).

    For NCBI, RefSeq and GenBank accession IDs are stored directly in this table.

    Updates the the genome dict "assembly_accession" field for each genome found.
    """
    # MySQL query
    database = "hgFixed"
    command = (
        "SELECT source,destination,matchCount "
        "FROM asmEquivalent "
        "WHERE sourceAuthority='ucsc' "
        "AND destinationAuthority!='ensembl' "
    )
    ret = query_ucsc(command, database)

    # convert to dataframe
    df = pd.DataFrame.from_records(ret)
    df.columns = ["name", "accession_name2", "match"]
    df.set_index("name", inplace=True)

    # get best match
    df["match_max"] = df.groupby(df.index)["match"].max()
    df = df[df["match"] == df["match_max"]]

    # extract accessions
    # example: GCF_000090745.1_AnoCar2.0 -> GCF_000090745.1
    df = df["accession_name2"].str.extract(r"(GC[AF]_\d{9}\.\d+)")

    # GCA > GCF
    duplicates = df.index.duplicated(False)
    unique_df = df[~duplicates]
    dup_df = df[duplicates]
    filtered_dup_df = dup_df[dup_df[0].str[2] == "A"]
    df = unique_df.append(filtered_dup_df)

    # During testing, 93/217 genomes we assigned an accession ID
    accession_series = df[0]
    for name, acc in accession_series.items():
        if name not in genomes:
            continue
        genomes[name]["assembly_accession"] = acc

    return genomes


def add_accessions2(genomes: dict) -> dict:
    """
    Some genomes have their assembly accession in the 'sourceName' field.

    Updates the the genome dict "assembly_accession" field for each genome found.
    """
    re_acc = re.compile(r"GC[AF]_\d{9}\.\d+")
    for name in genomes:
        hit = re_acc.search(genomes[name]["sourceName"])
        if hit:
            genomes[name]["assembly_accession"] = hit.group()
    return genomes


def add_annotation_links(genomes):
    """
    identify the available annotation types for each genome.
    """
    # MySQL query
    command = (
        "SELECT TABLE_SCHEMA,TABLE_NAME "
        "FROM information_schema.tables "
        "WHERE table_name IN "
        "('ensGene', 'ncbiRefSeq', 'knownGene', 'refGene') "
    )
    ret = query_ucsc(command, database=None)

    for name, annot in ret:
        if name not in genomes:
            continue
        genomes[name]["annotations"].append(annot)

    return genomes


def query_ucsc(command: str, database: str = None) -> Generator:
    """
    Execute a single MySQL query on the UCSC database.
    Streams the output into a generator.
    """
    cnx = mysql.connector.connect(
        host="genome-mysql.soe.ucsc.edu",
        user="genome",
        port=3306,
        database=database,
    )
    try:
        cur = cnx.cursor(buffered=False, raw=False)
        cur.execute(command)
        while True:
            ret = cur.fetchone()
            if not ret:
                break
            yield ret
    finally:
        cnx.close()


def download_annotation(name, annot, genomes_dir, localname, n=None):
    """
    Download the extended genePred file from the UCSC MySQL database.
    Next convert this to a BED and GTF file.
    """
    out_dir = os.path.join(genomes_dir, localname)
    mkdir_p(out_dir)
    tmp_dir = mkdtemp(dir=out_dir)
    pred_file = f"{os.path.join(tmp_dir, localname)}.annotation.extended.gp"
    gtf_file = f"{os.path.join(out_dir, localname)}.annotation.gtf"
    bed_file = f"{os.path.join(out_dir, localname)}.annotation.bed"

    # MySQL query 1: get column names for this genePred
    command = f"SHOW COLUMNS FROM {annot};"
    cols = list(query_ucsc(command, database=name))

    # drop columns the UCSC tools cannot handle
    # see https://genome.ucsc.edu/FAQ/FAQformat.html#format9
    accepted_cols = [
        "geneName",
        "name",
        "chrom",
        "strand",
        "txStart",
        "txEnd",
        "cdsStart",
        "cdsEnd",
        "exonCount",
        "exonStarts",
        "exonEnds",
        "score",
        "name2",
        "cdsStartStat",
        "cdsEndStat",
        "exonFrames",
    ]
    cols = [c[0] for c in cols if c[0] in accepted_cols]
    cols = ",".join(cols)

    # MySQL query 2: download genePred
    command = f"SELECT {cols} FROM {annot};"
    if n:
        command = f"SELECT {cols} FROM {annot} LIMIT {n};"
    ret = query_ucsc(command, database=name)

    # clean up genePred
    df = pd.DataFrame.from_records(ret)
    for c in [8, 9, 14]:
        if c in df:
            df[c] = df[c].str.decode("utf-8")
    df.to_csv(pred_file, index=False, header=False, sep="\t")

    # convert genePred to GTF and BED
    cmd = "genePredToGtf -source=genomepy file {0} {1}"
    sp.check_call(cmd.format(pred_file, gtf_file), shell=True)
    cmd = "genePredToBed {0} {1}"
    sp.check_call(cmd.format(pred_file, bed_file), shell=True)
    rm_rf(tmp_dir)


@cache
def scrape_accession(htmlpath: str) -> str:
    """
    Attempt to scrape the assembly accession (`GCA_`/`GCF_`) from a genome's readme.html,
    or any linked NCBI assembly pages can also be scraped.

    Parameters
    ----------
    htmlpath: str
        path to the readme.tml on hgdownload.soe.ucsc.edu

    Returns
    ------
    str
        Assembly accession or 'na'
    """
    ucsc_url = f"https://hgdownload.soe.ucsc.edu/{htmlpath}"
    try:
        text = read_url(ucsc_url)
    except (UnicodeDecodeError, urllib.error.URLError):
        return "na"

    # example accessions: GCA_000004335.1 (ailMel1)
    # regex: GC[AF]_ = GCA_ or GCF_, \d = digit, \. = period
    accession_regex = re.compile(r"GC[AF]_\d{9}\.\d+")
    match = accession_regex.search(text)
    if match:
        return match.group(0)

    # Search for an assembly link at NCBI
    match = re.search(r"https?://www.ncbi.nlm.nih.gov/assembly/\d+", text)
    if match:
        ncbi_url = match.group(0)
        try:
            text = read_url(ncbi_url)
        except (UnicodeDecodeError, urllib.error.URLError):
            return "na"

        # retrieve valid assembly accessions.
        # contains additional info, such as '(latest)' or '(suppressed)'. Unused for now.
        valid_accessions = re.findall(r"assembly accession:.*?GC[AF]_.*?<", text)
        text = " ".join(valid_accessions)
        match = accession_regex.search(text)
        if match:
            return match.group(0)
