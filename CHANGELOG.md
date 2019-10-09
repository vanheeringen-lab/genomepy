# Changelog

Here, the changes to `genomepy` will be summarized.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added

### Changed

### Fixed

## [0.6.0] - 2019-09-11

### Added

- Support for storing bzgip-compressed genomes (#41).

### Changed

- Removed support for Python 2 ([2020 is close!](https://python3statement.org/)).

### Fixed

- Ensembl annotation for non-vertebrate genomes should work again.
- Fixed bug where a deleted or empty config file would result in an error.

## [0.5.5] - 2019-03-19

### Added 

- Plugin for downloading genome blacklists (from Kundaje lab).

### Fixed

- Fix for new Ensembl REST API and FTP layout.
- Genomes from Ensembl with a space in their name can be downloaded.
- Plugin imports use relative parts to prevent conflicts with other imports.

## [0.5.4] - 2019-03-19

### Added

- Downloading annotation from NCBI now implemented.
- Genbank assemblies at NCBI can be searched and downloaded

### Fixed

- Fixed #23.
- Fixed #26.
- Fixed Ensembl downloads (#30)
- Fixed FTP tests for CI

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
