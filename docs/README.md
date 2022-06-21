:orphan:

# Generating the docs

## Run
* `sphinx-build docs build`
* open/refresh `build/index.html` in your browser (site did not refresh if the url ends with "#")
* if changes were made to the README.md in the root directory,
check that the lines still match in the various content pages.

## Sphinx setup
### install dependencies
* `mamba env update --file environment.yml`
* `mamba env update --file docs/requirements.yaml`
* `pip install -e .`

### follow guides:
* https://eikonomega.medium.com/getting-started-with-sphinx-autodoc-part-1-2cebbbca5365
(I manually moved all created files to ./docs)
* https://stackoverflow.com/a/62613202
(only change the marked lines in the template examples)

### changes to conf.py
* line 12: uncomment path setup
* line 15: `sys.path.insert(0, os.path.abspath('..'))`
* line 30: added extensions & extention configurations
* line 56: added `html_sidebars`
* line 65: `html_theme = "sphinx_rtd_theme"`
* line 66: added last-updated time
* line 73: added `html_context`

### entry point: `index.rst`
create the welcome page & Table of Content

### populate docs
add files to `content`
add figures to `images`

## Github-pages setup
* added the docs.yml to .github/workflows
* pushed to branch that triggers the workflow
* go to github > settings > pages
  * get a page and theme (required, but ignored by `.nojekyll`)
  * wait
