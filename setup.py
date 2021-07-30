import os
import sys

from setuptools import setup

# https://packaging.python.org/single_source_version/
exec(open("genomepy/__about__.py").read())

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()

with open("README.md", encoding="utf-8") as f:
    long_description = f.read()

# List of files and directories to include when packaging for release
packages = [
    "genomepy",
    "genomepy.annotation",
    "genomepy.config",
    "genomepy.genome",
    "genomepy.plugins",
    "genomepy.providers",
]
package_data = {"genomepy": ["config/default.yaml"]}
data_files = [("", ["LICENSE", "README.md", "CHANGELOG.md"])]

entry_points = {"console_scripts": ["genomepy=genomepy.cli:cli"]}

requires = [
    "appdirs",
    "biopython>=1.73",
    "bucketcache",
    "click",
    "colorama",
    "joblib",
    "loguru",
    "mygene",
    "mysql-connector-python",
    "norns>=0.1.5",
    "numpy",
    "pandas",
    "pyfaidx>=0.5.7",
    "requests",
    "tqdm>=4.51",
]

classifiers = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Developers",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Operating System :: MacOS :: MacOS X",
    "Operating System :: POSIX :: Linux",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name="genomepy",
    version=__version__,  # noqa: F821
    description="Automatic downloading and processing of genomes and metadata in command line and Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=packages,
    package_data=package_data,
    data_files=data_files,
    entry_points=entry_points,
    install_requires=requires,
    author=__author__,  # noqa: F821
    author_email="simon.vanheeringen@gmail.com",
    url="https://github.com/vanheeringen-lab/genomepy",
    license="MIT",
    classifiers=classifiers,
)
