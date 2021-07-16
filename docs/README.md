# Generating the docs

## Run
* install genomepy with `python setup.py develop` (needed only once)
* `sphinx-build docs build`
* open `build/index.html` in your browser

## Setup
### Install dependencies
* `mamba env update --file environment.yml`
* `mamba env update --file docs/requirements.yaml`

### follow steps on:
https://eikonomega.medium.com/getting-started-with-sphinx-autodoc-part-1-2cebbbca5365
(I manually moved all created files to ./docs)

https://stackoverflow.com/a/62613202
(only change the marked lines in the template examples)

### changes to conf.py
* line 12: uncomment path setup
* line 15: `sys.path.insert(0, os.path.abspath('..'))`
* line 30: added extensions & source_suffix
* line 54: `html_theme = "sphinx_rtd_theme"`
