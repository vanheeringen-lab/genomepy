# run:
$python setup.py develop
$cd .docs
$make html


# setup:
## dependencies
$mamba install sphinx sphinx_rtd_theme -y

## follow steps on:
https://eikonomega.medium.com/getting-started-with-sphinx-autodoc-part-1-2cebbbca5365
(I manually moved all created files to ./docs)

## changes to conf.py
line 15: # sys.path.insert(0, os.path.abspath('..'))
line 30: added 3 extensions
line 52: html_theme = "sphinx_rtd_theme"
