# dsRNAscan

[![CI Tests](https://github.com/rjandr/dsRNAscan/actions/workflows/ci-simple.yml/badge.svg)](https://github.com/rjandr/dsRNAscan/actions/workflows/ci-simple.yml)
[![Python](https://img.shields.io/badge/Python-3.8%2B%20(Linux)%20|%203.9%2B%20(macOS)-blue.svg)](https://www.python.org/downloads/)
[![Platforms](https://img.shields.io/badge/Platforms-Linux%20|%20macOS-green.svg)](https://github.com/Bass-Lab/dsRNAscan)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**dsRNAscan** is a bioinformatics tool for genome-wide identification of **double-stranded RNA (dsRNA) structures**. It uses a sliding window approach to detect inverted repeats that can form dsRNA secondary structures, with support for **G-U wobble base pairing** and **ML-based scoring**.

Browse human genome results at [dsrna.chpc.utah.edu](http://dsrna.chpc.utah.edu)

## Install

```bash
pip install dsrnascan
```

**Platforms:** Linux (Python 3.8+), macOS (Python 3.9+). Windows not supported (use WSL).

**Dependencies** (auto-installed): biopython, numpy, pandas, ViennaRNA

If ViennaRNA fails via pip, install with conda: `conda install -c bioconda viennarna`

## Usage

```bash
# Basic scan (defaults: -w 10000 -s 500 --score 75 -c 4)
dsrnascan input.fasta

# Specific chromosome with 8 CPUs
dsrnascan genome.fasta --only_seq chr21 -c 8

# Scan a specific region
dsrnascan genome.fasta --only_seq chr21 --start 33455482 --end 33655482

# Faster scan with larger step size (less overlap between windows)
dsrnascan genome.fasta -s 5000 -c 16

# Sensitive scan for shorter dsRNAs
dsrnascan sequence.fasta -w 5000 --min_bp 15 --paired_cutoff 60
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-w` | 10000 | Window size (bp) |
| `-s/--step` | 500 | Step size between windows |
| `-c/--cpus` | 4 | Number of CPUs |
| `--min_bp` | 25 | Minimum base pairs required |
| `--score` | 75 | Minimum einverted score |
| `--paired_cutoff` | 70 | Minimum % paired bases |
| `--only_seq` | None | Specific chromosome(s) to scan |
| `--start/--end` | Full seq | Region coordinates |
| `--forward-only` | False | Forward strand only |
| `--reverse-only` | False | Reverse strand only |
| `--no-ml` | False | Disable ML scoring |
| `--format` | bedpe | Output format: `bedpe`, `gff3`, or `both` |
| `--output-dir` | Auto | Output directory |

Run `dsrnascan --help` for advanced options (scoring parameters, folding temperature, etc.).

## Output

Results are written to the output directory:

**`*_merged_results.txt`** - Tab-delimited predictions with columns:

- **Coordinates:** Chromosome, Strand, i_start, i_end, j_start, j_end
- **einverted:** Score, RawMatch, PercMatch, Gaps
- **RNAduplex:** dG(kcal/mol), base_pairs, percent_paired, longest_helix, eff_i/j_start/end, i_seq, j_seq, structure
- **ML scores:** stability_model_score, probing_model_score, likely_edited (Yes/No), likely_forms (Yes/No)

**`*.bp`** - IGV arc visualization file

**`*.bedpe`** - BEDPE format with paired coordinates (one line per dsRNA) - default output

**`*.gff3`** - GFF3 with mRNA/exon types for IGV arc visualization (`--format gff3` or `--format both`)

## Browse Results

```bash
# Interactive viewer with Forna RNA structure visualization
dsrna-browse results_directory/

# With RNA editing site annotations (BED or GFF3)
dsrna-browse results_directory/ --editing-file editing_sites.bed
```

## Citation

If you use dsRNAscan, please cite:

> Comprehensive mapping of human dsRNAome reveals conservation, neuronal enrichment, and intermolecular interactions
> https://doi.org/10.1101/2025.01.24.634786

## Additional Tools

**overlap_analyzer** - Statistical enrichment analysis for genomic features overlapping dsRNA predictions. See [overlap_analyzer/README.md](overlap_analyzer/README.md). Not included in PyPI package; clone the repo to access.

## License

GNU General Public License v3.0 - see [LICENSE](LICENSE).

**Issues:** [GitHub Issues](https://github.com/Bass-Lab/dsRNAscan/issues)

## Acknowledgments

- EMBOSS team for the einverted algorithm
- ViennaRNA team for RNA folding algorithms
