import os.path
import re
import sys
from typing import Optional, Tuple, Iterable

from bisect import bisect
from glob import glob
from pyfaidx import Fasta, Sequence
from random import random
import pandas as pd
import numpy as np
from loguru import logger
import mygene

from genomepy.plugin import get_active_plugins
from genomepy.provider import ProviderBase
from genomepy.utils import (
    read_readme,
    write_readme,
    get_genomes_dir,
    get_localname,
    glob_ext_files,
    generate_fa_sizes,
    generate_gap_bed,
    safe,
)


class Genome(Fasta):
    """
    Get pyfaidx Fasta object of genome

    Also generates an index file of the genome

    Parameters
    ----------
    name : str
        Genome name

    genomes_dir : str
        Genome installation directory

    Returns
    -------
    pyfaidx.Fasta object
    """

    def __init__(self, name, genomes_dir=None):
        self.genomes_dir = get_genomes_dir(genomes_dir, check_exist=False)
        self.name = self._parse_name(name)
        self.filename = self._parse_filename(name)
        super(Genome, self).__init__(self.filename)

        # file paths
        self.genome_file = self.filename
        self.genome_dir = os.path.dirname(self.filename)
        self.index_file = self.filename + ".fai"
        self.sizes_file = self.genome_file + ".sizes"
        self.gaps_file = os.path.join(self.genome_dir, self.name + ".gaps.bed")
        self.readme_file = os.path.join(self.genome_dir, "README.txt")

        # genome attributes
        self.sizes = {}
        self.gaps = {}
        metadata = self._read_metadata()
        self.provider = metadata.get("provider")
        self.tax_id = metadata.get("tax_id")
        self.assembly_accession = metadata.get("assembly_accession")

    @property
    def sizes_file(self):
        return self.__sizes_file

    @sizes_file.setter
    def sizes_file(self, fname):
        """generate the sizes_file when the class is initiated"""
        if not os.path.exists(fname):
            generate_fa_sizes(self.genome_file, fname)
        self.__sizes_file = fname

    @property
    def sizes(self):
        return self.__sizes

    @sizes.setter
    def sizes(self, contig_sizes):
        self.__sizes = contig_sizes

    @sizes.getter
    def sizes(self):
        """Return sizes per contig when requested

        Returns
        -------
        contig_sizes : dict
            a dictionary with contigs as key and their lengths as values
        """
        if self.__sizes == {}:
            with open(self.sizes_file) as f:
                for line in f:
                    contig, length = line.strip().split("\t")
                    self.__sizes[contig] = int(length)
        return self.__sizes

    @property
    def gaps_file(self):
        return self.__gaps_file

    @gaps_file.setter
    def gaps_file(self, fname):
        """generate the gaps_file when the class is initiated"""
        if not os.path.exists(fname):
            generate_gap_bed(self.genome_file, fname)
        self.__gaps_file = fname

    @property
    def gaps(self):
        return self.__gaps

    @gaps.setter
    def gaps(self, gap_sizes):
        self.__gaps = gap_sizes

    @gaps.getter
    def gaps(self):
        """Return gap sizes per chromosome when requested

        Returns
        -------
        gap_sizes : dict
            a dictionary with chromosomes as key and the total number of
            Ns as values
        """
        if self.__gaps == {}:
            with open(self.gaps_file) as f:
                for line in f:
                    chrom, start, end = line.strip().split("\t")
                    start, end = int(start), int(end)
                    self.__gaps[chrom] = self.__gaps.get(chrom, 0) + end - start
        return self.__gaps

    @property
    def plugin(self):
        """dict of all active plugins and their properties"""
        p = dict()
        for plugin in get_active_plugins():
            p[plugin.name()] = plugin.get_properties(self)
        return p

    @property
    def annotation_gtf_file(self):
        return self.check_annotation_file("gtf")

    @property
    def annotation_bed_file(self):
        return self.check_annotation_file("bed")

    def _query_mygene(self, query: Iterable[str], fields: str = "genomic_pos"):
        # mygene.info only queries the most recent version of the Ensembl database
        # We can only safely continue if the local genome matched the Ensembl genome.
        # Even if the local genome was installed via Ensembl, we still need to check
        # if it is the same version
        result = self.ensembl_genome_info()
        if result is None:
            return None

        # Run the actual query
        logger.info("Querying mygene.info...")
        mg = mygene.MyGeneInfo()
        result = mg.querymany(
            query,
            scopes="symbol,name,ensembl.gene,entrezgene,ensembl.transcript,ensembl",
            fields=fields,
            species=self.tax_id,
            as_dataframe=True,
            verbose=False,
        )
        logger.info("Done")
        if "notfound" in result and result.shape[1] == 1:
            logger.error("No matching genes found")
            sys.exit()

        return result

    def ensembl_genome_info(self) -> Tuple[str, str, str]:
        """Return Ensembl genome information this genome.

        Returns
        -------
        (str, str, str)
            Ensembl name, accession, taxonomy_id
        """
        if self.provider == "Ensembl":
            return self.name, self.assembly_accession, self.tax_id

        # Fast lookup for some common queries.
        # We specifically provide the GRC genomes for mouse and human, as the patch
        # level does not influence the genome coordinates. However, they do have a
        # different assembly accession, so if we search by accession we don't get a
        # match.
        common_names = {
            "danRer11": "GRCz11",
            "hg38": "GRCh38.p13",
            "mm10": "GRCm38.p6",
            "dm6": "BDGP6.28",
        }
        if self.name in common_names:
            search_term = common_names[self.name]
        else:
            search_term = self.assembly_accession
            if search_term in ["na", None]:
                logger.warning(
                    "Cannot find a matching genome without an assembly accession."
                )
                return

        # search Ensembl by asssmbly accession or by specific Ensembl name (if we know it)
        logger.info(f"searching for {search_term}")
        results = list(ProviderBase.search_all(str(search_term), provider="Ensembl"))
        for name, _provider, accession, _species, tax_id, *_ in results:
            # Check if the assembly_id of the current Ensembl genome is the same as the
            # local genome. If it is identical, we can correctly assume that the genomes
            # sequences are identical.
            # For the genomes in the lookup table, we already know they match.
            if (
                common_names.get(self.name) == name
                or accession == self.assembly_accession
            ):
                return name, accession, tax_id

        logger.warning(
            "Could not find a matching assembly version in the current release of Ensembl."
        )

    def map_gene_dataframe(
        self, df: pd.DataFrame, genome: str, gene_field: str, product: str = "protein"
    ) -> pd.DataFrame:
        """Use mygene.info to map identifiers

        If the identifier can't be mapped, it will be dropped from the resulting
        annotation.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame with gene information. Should contain at least a "name" column.

        genome : str
            Genome name for mygene.info.

        gene_field : str
            Identifier for gene annotation. Uses mygene.info to map ids. Valid fields
            are: ensembl.gene, entrezgene, symbol, name, refseq, entrezgene. Note that
            refseq will return the protein refseq_id by default, use `product="rna"` to
            return the RNA refseq_id. Currently, mapping to Ensembl transcript ids is
            not supported.

        product : str, optional
            Either "protein" or "rna". Only used when `gene_field="refseq"`

        Returns
        -------
        pandas.DataFrame with mapped gene annotation.
        """
        cols = df.columns
        df["split_id"] = df["name"].str.split(r"\.", expand=True)[0]
        genes = df["split_id"].tolist()
        result = self._query_mygene(genes, fields=gene_field)
        if result is None:
            logger.error("Could not map using mygene.info")
        df = df.join(result, on="split_id")

        # Only in case of RefSeq annotation the product needs to be specified.
        if gene_field == "refseq":
            gene_field = f"{gene_field}.translation.{product}"

        # Get rid of extra columns from query
        df["name"] = df[gene_field]
        df = df[cols].dropna()

        return df

    def gene_annotation(
        self, gene_field: Optional[str] = None, product: str = "protein"
    ) -> pd.DataFrame:
        """Retrieve gene location from local annotation.

        You can use mygene.info to map identifiers on the fly by specifying
        `gene_field`. If the identifier can't be mapped, it will be dropped
        from the resulting annotation.

        Returns a DataFrame with gene annotation in bed12 format.

        Parameters
        ----------
        genome : str
            Genome name

        gene_field : str, optional
            Identifier for gene annotation. Uses mygene.info to map ids. Valid fields
            are: ensembl.gene, entrezgene, symbol, name, refseq, entrezgene. Note that
            refseq will return the protein refseq_id by default, use `product="rna"` to
            return the RNA refseq_id. Currently, mapping to Ensembl transcript ids is
            not supported.

        product : str, optional
            Either "protein" or "rna". Only used when `gene_field="refseq"`

        Returns
        -------
        pandas.DataFrame with gene annotation.
        """
        product = product.lower()
        if gene_field is not None:
            gene_field = gene_field.lower()

        if product not in ["rna", "protein"]:
            raise ValueError("Argument product should be either 'rna' or 'protein'")

        bed12_fields = [
            "chrom",
            "start",
            "end",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRrgb",
            "blockCount",
            "blockSizes",
            "blockStarts",
        ]

        if gene_field is not None:
            if gene_field == "refseq":
                anno_file = f"{self.name}.{gene_field}.{product}.annotation.bed"
            else:
                anno_file = f"{self.name}.{gene_field}.annotation.bed"

            target_anno = os.path.join(os.path.dirname(self.filename), anno_file)
            if os.path.exists(target_anno):
                return pd.read_csv(
                    target_anno,
                    sep="\t",
                    names=bed12_fields,
                    dtype={"chrom": "string", "start": np.uint32, "end": np.uint32},
                )

        bed = self.annotation_bed_file

        if bed is None:
            logger.info(f"No annotation file found for genome {self.name}!")
            logger.info("Run the following command to install the annotation:")
            logger.info(f"  genomepy install {self.name} {self.provider} --annotation")
            logger.info(
                f"Alternatively, copy your own annotation to {os.path.join(self.genome_dir, self.name + '.annotation.bed')}"
            )
            return

        df = pd.read_csv(bed, sep="\t", names=bed12_fields)

        # Optionally use mygene.info to map gene/transcript ids
        if gene_field is not None:
            logger.info("Mapping gene identifiers using mygene.info")
            logger.info(
                "This can take a long time, but results will be saved for faster access next time!"
            )
            df = self.map_gene_dataframe(
                df, self.name, gene_field=gene_field, product=product
            )

            # Save mapping results for quicker access next time
            df.to_csv(target_anno, sep="\t", index=False, header=False)

        return df

    @property
    def assembly_report(self):
        fname = os.path.join(self.genome_dir, "assembly_report.txt")
        if not os.path.exists(fname):
            if self.assembly_accession in ["na", None]:
                logger.warning(
                    "Can't download an assembly report without an assembly accession."
                )
                return
            logger.info("Assembly report not present, downloading...")
            p = ProviderBase.create("NCBI")
            p.download_assembly_report(self.assembly_accession, fname)

        return pd.read_csv(fname, sep="\t")

    @staticmethod
    def _parse_name(name):
        """extract a safe name from file path, url or regular names"""
        return os.path.basename(re.sub(".fa(.gz)?$", "", get_localname(name)))

    def _parse_filename(self, name):
        """
        accepts path to a fasta file, path to a fasta folder, or
        the name of a genome (e.g. hg38).

        returns the abspath to the fasta file
        """
        path_name = os.path.abspath(os.path.expanduser(name))
        if os.path.isfile(path_name):
            return path_name

        default_genome_dir = os.path.join(self.genomes_dir, self.name)
        for f in glob_ext_files(path_name) + glob_ext_files(default_genome_dir):
            if self.name + ".fa" in os.path.basename(f):
                return f

        raise FileNotFoundError(
            f"could not find {self.name}.fa(.gz) in genome_dir {default_genome_dir}"
        )

    def check_annotation_file(self, ext):
        """returns (gzipped) annotation file"""
        pattern = os.path.join(self.genome_dir, self.name + ".annotation.")
        file = glob(pattern + ext + "*")
        return file[0] if file else None

    @staticmethod
    def _update_provider(metadata):
        """check if provider is missing and try to update"""
        metadata["provider"] = "Unknown"
        url = metadata.get("genome url", "").lower()
        for provider in ["Ensembl", "UCSC", "NCBI"]:
            if provider.lower() in url:
                metadata["provider"] = provider
                break

    @staticmethod
    def _update_tax_id(metadata, provider=None, genome=None):
        """check if tax_id is missing and try to update"""
        taxid = "na"
        if genome:
            taxid = provider.genome_taxid(genome)
            taxid = str(taxid) if taxid != 0 else "na"
        metadata["tax_id"] = taxid

    @staticmethod
    def _update_assembly_accession(metadata, provider=None, genome=None):
        """check if assembly_accession is missing and try to update"""
        accession = "na"
        if genome:
            accession = provider.assembly_accession(genome)
        metadata["assembly_accession"] = accession

    def _update_metadata(self, metadata):
        """check if there is missing info that can be updated"""
        print("Updating metadata in README.txt", file=sys.stderr)
        if metadata.get("provider", "na") == "na":
            self._update_provider(metadata)

        known_provider = metadata["provider"] in ["Ensembl", "UCSC", "NCBI"]
        name = safe(metadata.get("original name", ""))
        missing_info = any(
            key not in metadata for key in ["tax_id", "assembly_accession"]
        )
        p = genome = None
        if known_provider and name and missing_info:
            p = ProviderBase.create(metadata["provider"])
            genome = p.genomes.get(name)

        if "tax_id" not in metadata:
            self._update_tax_id(metadata, p, genome)
        if "assembly_accession" not in metadata:
            self._update_assembly_accession(metadata, p, genome)

    def _read_metadata(self):
        """
        Read genome metadata from genome README.txt (if it exists).
        """
        metadata, lines = read_readme(self.readme_file)

        if (
            metadata.get("provider", "na") == "na"
            or "tax_id" not in metadata
            or "assembly_accession" not in metadata
        ) and os.access(self.readme_file, os.W_OK):
            self._update_metadata(metadata)
            write_readme(self.readme_file, metadata, lines)

        return metadata

    def _bed_to_seqs(self, track, stranded=False, extend_up=0, extend_down=0):
        bufsize = 10000
        with open(track) as fin:
            lines = fin.readlines(bufsize)
            for line in lines:
                if line.startswith("#") or line.startswith("track"):
                    continue

                vals = line.strip().split("\t")
                chrom, start, end = str(vals[0]), int(vals[1]), int(vals[2])
                name = f"{chrom}:{start}-{end}"

                # there might be more...
                starts = [start]
                ends = [end]

                # BED4: add name column to name
                if len(vals) >= 4:
                    name = " ".join((name, vals[3]))

                # BED5: check strandedness
                rc = False
                if stranded and len(vals) >= 6:
                    rc = vals[5] == "-"

                # BED12: get all blocks
                if len(vals) >= 12:
                    starts = [int(x) for x in vals[11].split(",")[:-1]]
                    sizes = [int(x) for x in vals[10].split(",")[:-1]]
                    starts = [start + x for x in starts]
                    ends = [start + size for start, size in zip(starts, sizes)]
                # convert to 1-based counting
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

                # load more lines if needed
                lines += fin.readlines(1)

    def _region_to_seq(self, region, extend_up=0, extend_down=0):
        chrom, coords = region.strip().split(":")
        start, end = [int(c) for c in coords.split("-")]
        start += 1
        start -= extend_up
        end += extend_down
        seq = self.get_seq(chrom, start, end)
        return seq.seq

    def _regions_to_seqs(self, track, extend_up=0, extend_down=0):
        if isinstance(track, list):
            for region in track:
                name = region.strip()
                seq = self._region_to_seq(name, extend_up, extend_down)
                yield Sequence(name, seq)
        else:
            with open(track) as fin:
                bufsize = 10000
                lines = fin.readlines(bufsize)
                for region in lines:
                    name = region.strip()
                    seq = self._region_to_seq(name, extend_up, extend_down)
                    yield Sequence(name, seq)

                    # load more lines if needed
                    lines += fin.readlines()

    @staticmethod
    def get_track_type(track):
        # region_p example: "chr1:123-456"
        region_p = re.compile(r"^([^\s]+):(\d+)-(\d+)$")
        if not isinstance(track, (str, bytes)) and isinstance(track, (list, tuple)):
            if isinstance(track[0], (str, bytes)) and region_p.search(track[0]):
                return "interval"
        with open(track) as fin:
            line = fin.readline().strip()
        if region_p.search(line):
            return "interval"
        return "bed"

    def track2fasta(
        self, track, fastafile=None, stranded=False, extend_up=0, extend_down=0
    ):
        """
        Return a list of fasta sequences as Sequence objects
        as directed from the track(s).

        Params:
        track: list/region file/bed file
            region(s) you wish to translate to fasta.
            Example input files can be found in genomepy/tests/data/regions.*

        fastafile: bool , optional
            return Sequences as list or save to file? (default: list)

        stranded: bool , optional
            return sequences for both strands? Required BED6 (or higher) as input (default: False)

        extend up/down: int , optional
            extend the sequences to either side? (command is strand sensitive, default: 0)
        """
        track_type = self.get_track_type(track)
        if track_type == "interval":
            seqqer = self._regions_to_seqs(
                track, extend_up=extend_up, extend_down=extend_down
            )
        else:
            seqqer = self._bed_to_seqs(
                track, stranded=stranded, extend_up=extend_up, extend_down=extend_down
            )

        if fastafile:
            with open(fastafile, "w") as fout:
                for seq in seqqer:
                    fout.write(f"{seq.__repr__()}\n")
        else:
            return [seq for seq in seqqer]

    @staticmethod
    def _weighted_selection(list_of_tuples, n):
        """
        Selects number random elements from a list of (weight, item) tuples.
        Based on code snippet by Nick Johnson
        """
        cuml = []
        items = []
        total_weight = 0.0
        for weight, item in list_of_tuples:
            total_weight += weight
            cuml.append(total_weight)
            items.append(item)

        return [items[bisect(cuml, random() * total_weight)] for _ in range(n)]

    def get_random_sequences(
        self, n=10, length=200, chroms=None, max_n=0.1, outtype="list"
    ):
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

        outtype : string , optional
            return the output as list or string.
            Options: "list" or "string", default: "list".

        Returns
        -------
        coords : list of lists/strings
            List with [chrom, start, end] genomic coordinates.
            String with "chrom:start-end" genomic coordinates
            (can be used as input for track2fasta).
        """
        if not chroms:
            chroms = self.keys()

        # dict of chromosome sizes after subtracting the number of Ns
        sizes = dict(
            [(chrom, len(self[chrom]) - self.gaps.get(chrom, 0)) for chrom in chroms]
        )

        # list of (tuples with) chromosomes and their size
        # (if that size is long enough for random sequence selection)
        lengths = [
            (sizes[x], x)
            for x in chroms
            if sizes[x] / len(self[x]) > 0.1 and sizes[x] > 10 * length
        ]
        if len(lengths) == 0:
            raise Exception("No contigs of sufficient size were found.")

        # random list of chromosomes from lengths (can have duplicates)
        chroms = self._weighted_selection(lengths, n)

        coords = []
        retries = 100
        cutoff = length * max_n
        for chrom in chroms:
            for _ in range(retries):
                start = int(random() * (sizes[chrom] - length))
                end = start + length
                count_n = self[chrom][start:end].seq.upper().count("N")
                if count_n <= cutoff:
                    break
            else:
                raise Exception(
                    f"Random subset ran {retries} times, "
                    f"but could not find a sequence with less than {cutoff} N's in {chrom}.\n"
                    "You can specify contigs using the CHROMS argument."
                )

            # list output example ["chr1", 123, 456]
            coords.append([chrom, start, end])

        if outtype != "list":
            # bed output example: "chr1:123-456"
            for i, region in enumerate(coords):
                coords[i] = [f"{region[0]}:{region[1]}-{region[2]}"]

        return coords
