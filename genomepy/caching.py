import os
from shutil import rmtree
from sqlite3 import DatabaseError
from time import sleep

from appdirs import user_cache_dir
from diskcache import Cache

from genomepy.__about__ import __version__

# Set short/long cache expiration times (in seconds)
cache_exp_short = 3.6e3
cache_exp_long = 8.64e4


genomepy_cache_dir = os.path.join(user_cache_dir("genomepy"), __version__)
os.makedirs(genomepy_cache_dir, exist_ok=True)

# Store the output of slow commands (marked with @disk_cache.memoize) for fast reuse
# DiskCache uses the LRU (least-recently-stored) eviction policy by default
try:
    disk_cache = Cache(directory=genomepy_cache_dir)
except DatabaseError:
    # another process was writing to the cache at the same time
    sleep(3)
    disk_cache = Cache(directory=genomepy_cache_dir)


def clean():
    """Remove cached data on providers."""
    rmtree(genomepy_cache_dir, ignore_errors=True)
    os.makedirs(genomepy_cache_dir, exist_ok=True)
    print("All clean!")
