import csv
import os
import subprocess as sp
from tempfile import mkdtemp

import numpy as np
import pandas as pd

from genomepy.files import gunzip_and_name, gzip_and_name
from genomepy.utils import rm_rf

GTF_FORMAT = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]

BED12_FORMAT = [
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


def _check_property(prop, fname):
    if prop is None:
        raise FileNotFoundError(f"'{fname}' required.")


def read_annot(fpath: str) -> pd.DataFrame:
    names = GTF_FORMAT if fpath.endswith(".gtf") else BED12_FORMAT
    df = pd.read_csv(
        fpath,
        sep="\t",
        names=names,
        dtype={"seqname": str, "start": np.uint32, "end": np.uint32},
        comment="#",
    )
    return df


def write_annot(df, fpath):
    df.to_csv(
        fpath,
        sep="\t",
        header=None,
        index=None,
        quoting=csv.QUOTE_NONE,
    )


def generate_annot(template, target, overwrite=False):
    """
    Create an annotation file type from the other file type.

    Parameters
    ----------
    template: str
        a GTF or BED filepath.

    target: str
        filepath to save the new annotation to.

    overwrite: bool, optional
        overwrite existing target file?
    """
    exts = os.path.basename(template.lower()).split(".")
    exts = [e for e in exts if e in ["gtf", "bed"]]
    if len(exts) == 0:
        raise ValueError("Template file must be in GTF or BED format.")
    template_ext = exts[-1]

    if not overwrite and os.path.exists(target):
        raise FileExistsError(f"{target} already exists! Set overwrite=True to ignore.")

    target_dir = os.path.dirname(target)
    tmp_dir = mkdtemp(dir=target_dir)
    tmp_target = os.path.join(tmp_dir, "new_annot")

    if template_ext == "bed":
        cmd = "bedToGenePred {0} /dev/stdout | genePredToGtf -source=genomepy file /dev/stdin {1}"
    else:
        cmd = "gtfToGenePred -ignoreGroupsWithoutExons {0} /dev/stdout | genePredToBed /dev/stdin {1}"

    # unzip template if needed
    template, is_unzipped = gunzip_and_name(template)
    # create new file
    sp.check_call(cmd.format(template, tmp_target), shell=True)
    # gzip if needed
    tmp_target = gzip_and_name(tmp_target, target.endswith(".gz"))
    _ = gzip_and_name(template, is_unzipped)

    os.replace(tmp_target, target)
    rm_rf(tmp_dir)


def bed_from_gtf(gtf_file, bed_file):
    """
    Create an annotation bed file from an annotation gtf file.
    Overwrites the existing bed file.
    """
    generate_annot(gtf_file, bed_file, overwrite=True)


def gtf_from_bed(bed_file, gtf_file):
    """
    Create an annotation gtf file from an annotation bed file.
    Overwrites the existing gtf file.
    """
    generate_annot(bed_file, gtf_file, overwrite=True)
