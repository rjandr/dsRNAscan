# Changelog

All notable changes to dsRNAscan will be documented in this file.

## [0.5.5] - 2026-07-31

### Fixed
- **Quality filters now run before (and again after) nested elimination.** Previously
  the order was dedup -> nested elimination -> filter, so a structure that would be
  filtered out could first eliminate a nested structure that passed the filters,
  losing both. Output is now a **strict superset** of the old behaviour - verified
  across 7 flag combinations, zero structures lost in any of them, all shared
  structures byte-identical in every column, recovering up to 6 structures per 100 kb.

### Changed
- **~22% faster with no loss of structures** (49.5s -> 38.7s on a 100 kb scan, 4 CPUs).
  Two changes, neither of which alters a surviving structure:
  - Hits below `--paired_cutoff` are now dropped **before** RNAduplex folding rather
    than after. `match_perc` comes from einverted, so it is known pre-fold; folding
    those structures and discarding them later was pure waste. They skew long-armed
    and `duplexfold` is O(n*m), so on a 100 kb test 20% of pairs consumed 69% of fold
    time.
  - Fold work is submitted longest-pair-first (LPT). Arm lengths vary ~100x, so
    submitting the largest first stops a late straggler from setting the makespan.
    `as_completed` reassembles by index, so ordering cannot affect output.
- Removed dead code: `predict_hybridization_batch` (defined, never called, and used a
  `ThreadPoolExecutor` that provides no parallelism because the ViennaRNA binding holds
  the GIL) and 6 unused imports.

### Added
- `--percent_paired_cutoff` separates the two metrics `--paired_cutoff` was silently
  gating. `--paired_cutoff` now applies only to einverted's match percentage; the new
  flag applies to RNAduplex's paired-base percentage and defaults to the same value,
  so default behaviour is unchanged. Set it to `0` to gate on the einverted match
  percentage only, as pre-0.4 releases did - on a test region that recovers ~6.6% more
  structures. Added for reproducing older result sets; it does not change any default.

### Changed
- **Faster einverted stage, identical results.** The bundled `einverted` now
  allocates its dynamic-programming tables once and reuses them across input
  records instead of re-allocating and O(maxrepeat^2)-zeroing them for every
  window. dsRNAscan groups each worker's windows into a single einverted call
  (multi-FASTA) to exploit this, so the DP matrix is initialised once per worker
  rather than once per window. Output is byte-identical to 0.5.4 (verified on
  ce11); parallelism and all result files are unchanged.
- **Fixed a memory leak exposed by the per-worker batching.** einverted's
  alignment-traceback scratch (`align1`/`align2`, ~`4*maxrepeat` ints) was not
  freed on the early exit taken when a repeat's traceback runs past the window
  edge. Harmless in 0.5.4 (one window per einverted process, reclaimed at exit),
  but with batching one process now scans hundreds of windows, so the leak
  accumulated to gigabytes and OOM-killed large scans. A single 300-window batch
  at the default band leaked ~2.5 GB (grew to 2.9 GB); now flat at 381 MB (the DP
  matrix alone). Output is byte-identical.
- einverted source is now vendored and self-contained in `einverted_src/`
  (no EMBOSS/libajax build required); the G-U wobble patch is baked in. The
  Linux x86_64 binary was rebuilt from it; other platform binaries should be
  rebuilt from the same source for the release (see `einverted_src/README.md`).
- **Bounded peak memory.** Results are now finalized (deduplicate +
  nested-elimination + filter) and flushed **per chromosome** instead of
  accumulating the whole genome in one DataFrame, so peak memory scales with the
  largest chromosome rather than the whole genome. This is output-identical
  because the dedup key already includes `chromosome` and nested-elimination
  already groups by `(chromosome, strand)`  -  nothing is cross-chromosome. On a
  9 Mb multi-chromosome test, peak RSS dropped 260 -> 200 MB; the reduction grows
  with genome size and removes the out-of-memory wall that blocked large
  vertebrate genomes from scanning in-core.
- Sequence/structure columns use PyArrow-backed strings and coordinates use
  `int32` when available (best-effort; falls back silently if PyArrow is
  absent), shrinking the retained results further. Requires no new hard dependency.

## [0.5.4] - 2026-07-03

### Fixed
- `--max_span` can now restrict the search span below the window size. It previously used `max(window_size, max_span)`, so any value below `-w` was silently discarded and the flag was a no-op for its main use case (short-duplex searches). Now uses `min(window_size, max_span)`. **Results from earlier runs that set `--max_span` below `-w` did not have that constraint applied.**
- `eliminate_nested_dsrnas` rewritten as a sorted sliding-window sweep. It previously allocated ~7 pairwise `(n, n)` boolean matrices per chromosome-strand group (250 GB per matrix at n = 500k) and OOM'd on dense scans. Now O(n) memory with no pairwise matrices; ~110x faster at n = 10k, and 200k structures process in ~1.2s.
- Nested elimination writeback no longer calls `index.get_loc()` per row (was O(n^2) on a non-unique index); uses an O(1) positional map.
- Nested elimination logs a warning on pathologically dense loci (max comparison window > 5000) suggesting `--no-eliminate-nested`.

### Notes
- `--start`/`--end` verified to produce correct genomic coordinates on both strands. Region scans are not bit-identical to full scans near region boundaries (fewer overlapping windows, and einverted arm boundaries can vary by ~1 nt with window framing). Pad regions by at least one window width and filter afterward.

## [0.5.3] - 2026-07-03

### Removed
- `--min_len` and `--max_len` arm length filters (introduced in 0.5.1). Post-hoc filtering saved no compute over filtering the output directly, and the default `--min_len` silently interacted with `--min_bp`. Filter the TSV/BEDPE output directly if arm length limits are needed.

### Fixed
- `--no-eliminate-nested` now actually disables nested dsRNA removal (flag value was never passed to worker processes via ProcessorArgs)
- `--score`-only mode now correctly disables base-pair filtering (the `--min_bp` default of 25 previously made the score-only code path unreachable, silently applying a 25 bp floor)

## [0.5.2] - 2026-07-03

### Fixed
- `--min_len` and `--max_len` now correctly filter by arm length post-hoc (einverted's `-maxrepeat` flag controls total span, not arm length, and has no `-minrepeat` equivalent)
- Removed spurious `-minrepeat` flag from einverted command (was silently ignored)

## [0.5.1] - 2026-07-02

### Changed
- BEDPE is now the default output format (`--format bedpe|gff3|both`)
- `--min` and `--max` renamed to `--min_len` and `--max_len` to avoid confusion with `--min_bp`
- `--max_len` now defaults to window size (previously hardcoded 10000)
- `--eliminate-nested` replaced by `--no-eliminate-nested` (flag was previously non-functional)
- dsRNA browser updated with cleaner interface

### Fixed
- `--min_len` / `--max_len` now correctly passed to einverted as `-minrepeat` / `-maxrepeat` (previously silently ignored)
- `--max_span` now correctly passed through to worker processes via ProcessorArgs (previously silently ignored)
- Removed `--algorithm` and `--clean` arguments that had no effect

## [0.5.0] - 2026-04-04

### Added
- GFF3 and BEDPE output files generated automatically
- `likely_edited` and `likely_forms` columns based on ML model Youden thresholds
- DBN copy button in dsrna-browse for quick structure export
- Confidence badges in dsrna-browse (green "Likely Edited", blue "Likely Forms")
- BP files split by strand (forward.bp and reverse.bp) for cleaner IGV visualization

### Changed
- Default step size increased from 150 to 500 (20x window overlap, 3x faster scans)
- Renamed ML columns: `stability_score` to `stability_model_score`, `probing_score` to `probing_model_score`
- Corrected stability model Youden threshold to 0.2471
- ydf and scikit-learn are now required dependencies
- Overlap analyzer reads dsRNAscan TSV output directly, vectorized BED creation

### Fixed
- Sequences shorter than window size now create one full-length window
- Removed hardcoded paths from overlap analyzer

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