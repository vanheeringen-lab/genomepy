# Changelog

Here, the changes to `genomepy` will be summarized.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

## [Unreleased]

## [0.8.1] - 2020-05-11

### Added
- Now using the UCSC rest API
- `genomepy search` now accepts taxonomy IDs
- `genomepy search` will now return taxonomy IDs and Accession numbers
- The README.txt will now store taxonomy IDs and Accession numbers
- Gene annotations:
    - Downloading of annotation file (BED/GTF/GFF3) from URL
    - Automatic search for annotation file (GTF/GFF3) in genome directory when downloading from URL
    - Option for URL provider to link to annotation file (to process is similarly to other providers)
    - Automatic annotation sanitizing (and skip sanitizing flag `-s` for `genomepy install`)
    - Option to only download annotation with `genomepy install -o`
- Plugins:
    - Blacklists are automatically unzipped.
    - Multithreading support for plugins, thanks to alienzj!
    - STAR now uses the annotation file for a (one-pass) splice-aware index

### Changed
- sizes no longer a plugin, but always gets executed
- `genomepy FUNCTION --help` texts expanded
- all genomepy classes exported when imported into Python
- all providers now let you know when they are downloading assembly information.
- more descriptive feedback to installing & many errors

### Removed
- Sizes plugin
- Old tests
- Removed outdated dependency `xmltodict`

### Fixed
- `genomepy config` options made more robust
- README.txt will no longer:
  - update 3x for each command
  - drop regex info
  - have duplicate lines

### Refactoring
- Genome class moved to `genome.py`
- Many functions moved to `utils.py`
- Many other functions made static methods of a class
- `Genome.track2fasta` and `Genome.get_random_sequence` optimized
- All Provider classes now store their genomes as a dict-in-dict, with the assembly name as key.
- Many Provider class functions now standardized. Many functions moved to from the daughter classes to the ProviderBase class.
- README.txt file generation and updating standardized
- Unit tests! all functions now have an individual test. Almost all test use functions already tested prior to them.
- Old tests incorporated in several extra tests (e01, e02, e03).
- Raise statements now use more fitting errors
- All instances of `os.remove` exchanged for `os.unlink`
- Almost all warnings fixed
- Extensive, COVID19-enabled, and somewhat pointless alphabetizing, optimizing and/or organizing changes to
    - imports everywhere
    - `.gitignore`
    - `.travis.yml`
    - `release_checklist.md`
    - `cli.py`
    - strings (many strings with .format() replaced with f-strings)

## [0.7.2] - 2019-03-31

### Fixed
- Fix minor issue with hg19 wrong blacklist url
- Ensembl downloads over http instead of https (release 99 no longer has https)

## [0.7.1] - 2019-11-20

### Fixed
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

- Fixed genomes_dir argument to `genomepy install`
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
