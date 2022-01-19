import os
from shutil import rmtree

from appdirs import user_cache_dir
from diskcache import Cache

from genomepy.__about__ import __version__

# Cache expiration times (in seconds)
cache_exp_short = 3600
cache_exp_long = 3600 * 24

genomepy_cache_dir = os.path.join(user_cache_dir("genomepy"), __version__)
os.makedirs(genomepy_cache_dir, exist_ok=True)

# Store the output of slow commands (marked with @disk_cache.memoize) for fast reuse.
# diskcache uses the LRU eviction policy by default
disk_cache = Cache(directory=genomepy_cache_dir)


def clean():
    """Remove cached data on providers."""
    rmtree(genomepy_cache_dir, ignore_errors=True)
    os.makedirs(genomepy_cache_dir, exist_ok=True)
    print("All clean!")
