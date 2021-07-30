import os
from shutil import rmtree

from appdirs import user_cache_dir
from bucketcache import Bucket
from joblib import Memory

from genomepy.__about__ import __version__

genomepy_cache_dir = os.path.join(user_cache_dir("genomepy"), __version__)
os.makedirs(genomepy_cache_dir, exist_ok=True)

# Store the output of slow commands (marked with @cache and @goldfish_cache) for fast reuse.
# Bucketcache creates a new pickle for each function + set of unique variables,
# For class methods, use ignore=["self"] or use @staticmethod.
cache = Bucket(genomepy_cache_dir, days=7)
goldfish_cache = Bucket(genomepy_cache_dir, minutes=10)

memory = Memory(genomepy_cache_dir, verbose=0)


def clean():
    """Remove cached data on providers."""
    rmtree(genomepy_cache_dir, ignore_errors=True)
    os.makedirs(genomepy_cache_dir, exist_ok=True)
    print("All clean!")
