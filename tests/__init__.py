"""
Global variables for pytest
"""
import os
from platform import system

linux = system() == "Linux"
travis = os.environ.get("TRAVIS") == "true"
