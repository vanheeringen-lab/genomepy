import os
import sys
import pypandoc

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

# https://packaging.python.org/single_source_version/
exec(open("genomepy/__about__.py").read())

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()

# trick to get rst file for PyPi, see:
# http://stackoverflow.com/questions/26737222/pypi-description-markdown-doesnt-work/26737672#26737672
# For the recent versions of pandoc, pypandoc >=1.3.3 is needed.
long_description = pypandoc.convert("README.md", "rst")

packages = ["genomepy", "genomepy/plugins"]

package_data = {"genomepy": ["cfg/*.yaml"]}

requires = [
    "pytest",
    "click",
    "pyfaidx>=0.5.1",
    "norns>=0.1.5",
    "xmltodict",
    "bucketcache",
    "requests",
    "appdirs",
    "psutil",
]

entry_points = {"console_scripts": ["genomepy=genomepy.cli:cli"]}

classifiers = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Developers",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Operating System :: MacOS :: MacOS X",
    "Operating System :: MacOS :: MacOS X",
    "Operating System :: POSIX :: Linux",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.3",
    "Programming Language :: Python :: 3.4",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name="genomepy",
    version=__version__,  # noqa: F821
    description="Genomes in Python",
    long_description=long_description,
    packages=packages,
    package_data=package_data,
    entry_points=entry_points,
    install_requires=requires,
    author=__author__,  # noqa: F821
    author_email="simon.vanheeringen@gmail.com",
    url="https://github.com/simonvh/genomepy",
    license="MIT",
    classifiers=classifiers,
)
