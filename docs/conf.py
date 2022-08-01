# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

project = 'genomepy'
copyright = 'Siebren Frölich, Maarten van der Sande, Tilman Schäfers and Simon van Heeringen'
author = 'Siebren Frölich, Maarten van der Sande, Tilman Schäfers and Simon van Heeringen'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',      # automatic documentation from docstrings
    'sphinx.ext.coverage',     # gather documentation coverage stats
    'sphinx.ext.napoleon',     # recognize numpy & google style docstrings
    'sphinx.ext.autosummary',  # Create neat summary tables
    'm2r2',                    # recognize markdown files
]

# Configuration of sphinx.ext.autosummary
autosummary_generate = True

# Configuration of m2r2
source_suffix = ['.rst', '.md']

# Configuration of sphinx.ext.coverage
coverage_show_missing_items = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# add the global ToC to the sidebar of each page
html_sidebars = {
    "**": ["globaltoc.html", "relations.html", "sourcelink.html", "searchbox.html"]
}

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
html_last_updated_fmt = '%Y-%m-%d, %H:%M (UTC)'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# add the "Edit on GitHub" link to the docs
html_context = {
    "display_github": True,
    "github_host": "github.com",
    "github_user": "vanheeringen-lab",
    "github_repo": "genomepy",
    "github_version": "develop",
    "conf_py_path": "/docs/",
    "source_suffix": [".rst", ".md"],
}
