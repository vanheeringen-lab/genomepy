import os
import time
from shutil import rmtree
from functools import wraps

from appdirs import user_cache_dir
from bucketcache import Bucket
from joblib import Memory,dump,load,hash
from diskcache import Cache

from genomepy.__about__ import __version__

genomepy_cache_dir = os.path.join(user_cache_dir("genomepy"), __version__)
os.makedirs(genomepy_cache_dir, exist_ok=True)
# Store the output of slow commands (marked with @cache and @goldfish_cache) for fast reuse.
# Bucketcache creates a new pickle for each function + set of unique variables,
# For class methods, use ignore=["self"] or use @staticmethod.
cache = Bucket(genomepy_cache_dir, days=7)
goldfish_cache = Bucket(genomepy_cache_dir, minutes=10)

memory = Memory(genomepy_cache_dir, verbose=0)
disk_cache = Cache(directory=genomepy_cache_dir)


def clean():
    """Remove cached data on providers."""
    rmtree(genomepy_cache_dir, ignore_errors=True)
    os.makedirs(genomepy_cache_dir, exist_ok=True)
    print("All clean!")


def hybcache(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        # Define short/long term cache expiration in seconds
        exp_short = 40 # type: int
        exp_long = 3600 # type: int   
        # Generate the cache key from the function's arguments.
        kwd_mark = object()
        key = args + (kwd_mark,) + tuple(sorted(kwargs.items()))
        print(key)
        key = hash(key)
        # Pickle object
        ext = ".joblib"
        pickle_obj = os.path.join(disk_cache.directory,f"{key}{ext}") #Check if value exists in the short term cache
        value = disk_cache.get(key)
        #Check if the function call deos not exists in the memory cache
        if value is None:
            #Check if the value exists as pickle and if it is expired
            if os.path.exists(pickle_obj):
                current_time = time.time()
                expiry_time = os.path.getctime(pickle_obj) + exp_long
                # Remove pickled object (long-term) if expired 
                if (current_time > expiry_time):
                    print(f"Removing {key} from short-term memory since it is expired")
                    os.remove(pickle_obj)
                    # We need to re-compute and store the result in memory and as pickle
                    print(f"re-computing {key} and caching value")
                    value = func(*args, **kwargs)
                    print(value)
                    disk_cache.close()
                    with Cache(disk_cache.directory) as reference:
                        reference.set(key, value, exp_short)
                    # Store function result as pickle
                    dump(value, os.path.join(disk_cache.directory,pickle_obj), protocol=3, cache_size=2000, compress=True)
                else:
                    # Load existing pickle into memory and skip computation
                    print(f"Loading {key} from pickle into cache")
                    with open(pickle_obj, 'rb') as pickle_handle:
                        disk_cache.close() 
                        with Cache(disk_cache.directory) as reference:
                            reference.set(key, load(pickle_obj), exp_short)
                    value = disk_cache.get(key)
            else:
                # Run the function and cache the result for next time.
                print(f"computing {key} and storing value inner")
                value = func(*args, **kwargs)
                print(value)
                disk_cache.close()
                with Cache(disk_cache.directory) as reference:
                    reference.set(key, value, exp_short)
                # Store function result as pickle
                dump(value, os.path.join(disk_cache.directory,pickle_obj),  protocol=3, cache_size=2000, compress=True)
        else:
            print(f"Loading {key} from short-term memory")
            # Skip the function entirely and use the cached value instead.
            value= disk_cache.get(key)
        return value
    return wrapper