import pytest

import genomepy
from genomepy.caching import disk_cache, clean
from genomepy.annotation import query_mygene


def test_caching_query_mygene():
    # Expected caching key for GRCz11 annotation
    a = genomepy.Annotation("GRCz11", "data")
    caching_key = (
        "genomepy.annotation.mygene.query_mygene",
        a.genes(),
        7955,
        "symbol",
        None,
    )
    x = query_mygene(a.genes(), 7955, "symbol")
    # Retrieved cached results for query_mygene
    y = disk_cache.get(caching_key)
    x = x.rename_axis("genes").reset_index()
    y = y.rename_axis("genes").reset_index()
    # Delete key from cache object
    disk_cache.delete(caching_key)
    # Check that results before/after caching are identical
    assert (
        x.equals(y[x.columns]) == True
    ), "Cached query_mygene output does not match query output"
