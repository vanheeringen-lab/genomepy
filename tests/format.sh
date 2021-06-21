#!/usr/bin/env bash

: '
Command line script run formatters and linters.
Run from root directory with `  bash tests/format.sh  `

Formatters:
- isort
- autoflake
- black

Linters:
- flake8 + flake8-bugbear
'

targets=("setup.py genomepy/ tests/")

autoflake \
  --recursive \
  --remove-all-unused-imports \
  --remove-duplicate-keys \
  --remove-unused-variables \
  --ignore-init-module-imports \
  $targets

isort \
  --overwrite-in-place \
  --profile black \
  --conda-env environment.yml \
  $targets


black \
  --quiet \
  $targets


flake8 \
  $targets


git status
