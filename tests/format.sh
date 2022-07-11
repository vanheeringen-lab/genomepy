#!/usr/bin/env bash

: '
Format or lint the repo.

Run from root directory to format:
  bash tests/format.sh
Run from root directory to lint:
  bash tests/format.sh lint
'

# files & directories to format/lint
targets=("genomepy/ tests/")

# check current directory
if [[ $(pwd) != *genomepy ]] || [[ $(ls) != *genomepy* ]]; then
  echo "Script must be run from the base of the genomepy repo!"
  exit 1
fi

# parse arg
lint=false
if [ "$1" = "lint" ]; then
  lint=true
fi

# store diff before formatting
before=$(git diff)


autoflake \
  $( $lint && echo '--check' ) \
  --in-place \
  --recursive \
  --remove-all-unused-imports \
  --remove-duplicate-keys \
  --remove-unused-variables \
  $targets \
  | grep -v 'No issues detected!'

isort \
  $( $lint && echo '--check' ) \
  --overwrite-in-place \
  --profile black \
  --conda-env environment.yml \
  $targets

black \
  $( $lint && echo '--check' ) \
  $targets \
  2>&1 \
  | grep 'would reformat'

if $lint; then
  flakeheaven lint \
  $targets

  echo ""
  echo "Done"
  exit 0
fi


# show formatted files
after=$(git diff)
if [ "$before" = "$after" ]; then
  echo "No changes made!"
else
  echo ""
  git status
fi
