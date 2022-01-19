from time import sleep
import pytest

import genomepy
from genomepy.caching import disk_cache, clean
from genomepy.annotation import query_mygene


def test_caching_query_mygene():
    # Create test annotation
    a = genomepy.Annotation("GRCz11", genomes_dir="tests/data")
    # Expected caching key for GRCz11 annotation
    caching_key = (
        "genomepy.annotation.mygene.query_mygene",
        a.genes(),
        7955,
        "symbol",
        None,
    )
    # Store output for comparision with cached data
    x = query_mygene(a.genes(), 7955, "symbol")
    x = x.rename_axis("genes").reset_index()
    # Retrieved cached results for query_mygene
    y = disk_cache.get(caching_key)
    y = y.rename_axis("genes").reset_index()
    # Check that results before/after caching are identical
    assert (
        x.equals(y[x.columns]) == True
    ), "Cached query_mygene output does not match query output"
    # Delete key from cache object
    disk_cache.delete(caching_key)
