# Changelog

Here, the changes to `genomepy` will be summarized.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added

- Downloading annotation from NCBI now implemented.
- Genbank assemblies at NCBI can be searched and downloaded

### Fixed

- Fixed #23.
- Fixed #26.

## [0.5.2] - 2018-09-11

### Fixed

- Fixed genome_dir argument to `genomepy install`
- Fixed msgpack dependency
- Fixed issue with `config generate` where config directory does note exist.

## [0.3.1]

- Added requests dependency
- Removed dependency on xdg, as it didn't support OSX
- Fixed string decoding bug

## [0.3.0]

- Started CHANGELOG.
- Genome listings are cached locally.
- Added `-m hard` option to `install` to hard-mask sequences.
- Added `-l` option to `install` for a custom name.
- Added `-r` and `--match/--no-match` option to select sequences by regex.
