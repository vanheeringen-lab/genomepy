import os
from shutil import rmtree
from time import time

from appdirs import user_cache_dir
from diskcache import Cache
from filelock import FileLock

from genomepy.__about__ import __version__
from genomepy.utils import rm_rf

# Set short/long cache expiration times (in seconds)
cache_exp_short = 3.6e3
cache_exp_long = 8.64e4


genomepy_cache_dir = os.path.join(user_cache_dir("genomepy"), __version__)
os.makedirs(genomepy_cache_dir, exist_ok=True)

# create a lock, so only one tread can access the cache at once
lock_file = os.path.join(genomepy_cache_dir, "cache.lock")
if os.path.exists(lock_file) and time() - os.stat(lock_file).st_mtime > 60:
    # remove abandoned lock
    rm_rf(lock_file)
lock = FileLock(lock_file)

with lock:
    # Store the output of slow commands (marked with @disk_cache.memoize) for fast reuse
    # DiskCache uses the LRU (least-recently-stored) eviction policy by default
    disk_cache = Cache(directory=genomepy_cache_dir)


@lock
def clean():
    """Remove cached data on providers."""
    rmtree(genomepy_cache_dir, ignore_errors=True)
    os.makedirs(genomepy_cache_dir, exist_ok=True)
    print("All clean!")
