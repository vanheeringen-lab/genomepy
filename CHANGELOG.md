# Changelog

Here, the changes to `genomepy` will be summarized.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

## [Unreleased]

## [0.7.2] - 2019-03-31

### Fixes
- Fix minor issue with hg19 wrong blacklist url
- Ensembl downloads over http instead of https (release 99 no longer has https)

## [0.7.1] - 2019-11-20

### Fixes
- STAR is not longer enabled by default

## [0.7.0] - 2019-11-18

### Added
- Direct downloading from url through url provider.
- Added `--force` flag. Files will no longer be overwritten by default.
- Provider specific options:
  - `--ensembl-version`: specify release version.
  - `--ensembl-toplevel`: by default, `genomepy install` will search for primary assemblies. 
  This flag will only download toplevel assemblies.
- Added STAR index plugin

### Changed
- Providers are now case-insensitive.
- Extended testing.
- Increased minimal Python version to 3.6.
- Removed gaps from plugins, added gaps to core functionality.

### Fixed
- bugfix: NCBI will show all versions of an assembly (will no longer filter on BioSample ID, instead filters on asm_name).
- fix: gaps file will be generated when needed.

## [0.6.1] - 2019-10-10

### Fixed

- Fixed bug with get_track_type.

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
