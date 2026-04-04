# Changelog

All notable changes to dsRNAscan will be documented in this file.

## [0.4.8] - 2026-04-04

### Added
- GFF3 and BEDPE output formats generated automatically
- `likely_edited` and `likely_forms` categorical columns based on ML model thresholds
- ML models bundled in package (`dsrnascan/models/`)
- ydf and scikit-learn now required dependencies (no longer optional)
- Overlap analyzer supports dsRNAscan TSV output directly (not just parquet)

### Changed
- Renamed ML columns: `stability_score` to `stability_model_score`, `probing_score` to `probing_model_score`
- Corrected stability model Youden threshold (0.2471)
- Overlap analyzer: vectorized BED creation, removed hardcoded paths
- Removed `conda-recipe/`, `test_data/`, and einverted source files from repo

## [0.4.7] - 2026-04-04

### Added
- ML scoring with stability and probing models (disable with `--no-ml`)
- tqdm progress bars during einverted processing
- `dsrna-browse`: server-side pagination, sortable/filterable table, sequence copy buttons, arm length columns
- Forna rendering size limit with fallback message for large structures

### Changed
- Deduplicate einverted results before RNAduplex (~10x fewer calls)
- Removed legacy/batch code paths and duplicate logging
- Cleaned up README with concise parameter table

### Fixed
- `--start`/`--end` now correctly trims sequence before windowing
- No more truncated/partial windows at sequence boundaries
- Exit code issue where result count was returned as exit code
- SettingWithCopyWarning in RNAduplex deduplication

## [0.4.6] - 2025-09-15

### Changed
- **Major Update**: dsRNAscan no longer requires EMBOSS installation
- Includes pre-compiled einverted binaries for all major platforms
- Platform-specific binary selection: Linux (x86_64, ARM64, i386), macOS (x86_64, ARM64), Windows (x86_64)
- Simplified installation to just `pip install dsrnascan`

### Added
- Standalone einverted binary that doesn't require EMBOSS libraries
- G-U wobble base pairing support in einverted
- Platform-specific wheel distribution on PyPI
- Automatic platform detection and binary selection

### Fixed
- Score header output compatibility with dsRNAscan parser
- stdin/stdout communication between dsRNAscan and einverted
- Binary permissions automatically set during installation

### Removed
- Dependency on EMBOSS installation
- Complex compilation steps from installation process
- Generic einverted fallback in favor of platform-specific binaries

## [0.4.5] - 2025-09-01

### Changed
- Improved DataFrame processing performance
- Enhanced memory efficiency for large-scale analyses

## [0.4.4] - 2025-08-15

### Added
- Initial support for chunked processing
- Improved progress reporting

### Fixed
- Memory issues with very large FASTA files