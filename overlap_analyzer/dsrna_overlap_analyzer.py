#!/usr/bin/env python3
"""
Universal dsRNA Overlap Analyzer
Modular tool for analyzing overlap between any genomic dataset and dsRNA regions

Usage:
    python dsrna_overlap_analyzer.py input.bed [options]
"""

import argparse
import sys
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import pandas as pd
import numpy as np
import pybedtools
from abc import ABC, abstractmethod
import json
import logging
from datetime import datetime
from tqdm import tqdm
import tempfile
import atexit
import shutil

# Try to import enhanced output module
try:
    from enhanced_output import (
        print_enhanced_results, 
        create_enrichment_plot,
        print_subset_recommendation
    )
    HAS_ENHANCED_OUTPUT = True
except ImportError:
    HAS_ENHANCED_OUTPUT = False

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Global temp directory for this session
TEMP_DIR = None

def setup_temp_dir():
    """Setup a session-specific temp directory"""
    global TEMP_DIR
    TEMP_DIR = tempfile.mkdtemp(prefix="dsrna_analyzer_")
    logger.debug(f"Created temp directory: {TEMP_DIR}")
    
    # Register cleanup function
    atexit.register(cleanup_temp_dir)
    
def cleanup_temp_dir():
    """Clean up all temp files at exit"""
    global TEMP_DIR
    if TEMP_DIR and os.path.exists(TEMP_DIR):
        shutil.rmtree(TEMP_DIR)
        logger.debug(f"Cleaned up temp directory: {TEMP_DIR}")

def get_temp_filename(prefix="temp", suffix=".bed"):
    """Get a temp filename in our session directory"""
    if not TEMP_DIR:
        setup_temp_dir()
    return os.path.join(TEMP_DIR, f"{prefix}_{os.getpid()}_{np.random.randint(100000)}{suffix}")


class InputParser(ABC):
    """Abstract base class for input file parsers"""
    
    @abstractmethod
    def parse(self, filepath: str) -> pybedtools.BedTool:
        """Parse input file and return BedTool object"""
        pass
    
    @abstractmethod
    def get_metadata(self) -> Dict:
        """Extract metadata from input file"""
        pass


class BEDParser(InputParser):
    """Parser for BED format files"""
    
    def parse(self, filepath: str) -> pybedtools.BedTool:
        """Parse BED file"""
        return pybedtools.BedTool(filepath)
    
    def get_metadata(self) -> Dict:
        """Extract metadata from BED file"""
        return {
            'format': 'BED',
            'strand_aware': True  # Assumes column 6 has strand info
        }


class GFF3Parser(InputParser):
    """Parser for GFF3 format files"""
    
    def __init__(self, feature_types: Optional[List[str]] = None):
        self.feature_types = feature_types
    
    def parse(self, filepath: str) -> pybedtools.BedTool:
        """Parse GFF3 file and convert to BED format"""
        features = []
        
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                chrom, source, feature_type, start, end, score, strand, phase, attributes = parts
                
                # Filter by feature type if specified
                if self.feature_types and feature_type not in self.feature_types:
                    continue
                
                # Convert to 0-based coordinates
                start = int(start) - 1
                end = int(end)
                
                # Extract ID from attributes
                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value
                
                feature_id = attr_dict.get('ID', f"{feature_type}_{chrom}:{start}-{end}")
                
                features.append([chrom, start, end, feature_id, score, strand])
        
        # Create temporary BED file
        temp_bed = get_temp_filename("gff3", ".bed")
        pd.DataFrame(features).to_csv(temp_bed, sep='\t', header=False, index=False)
        
        return pybedtools.BedTool(temp_bed)
    
    def get_metadata(self) -> Dict:
        """Extract metadata from GFF3 file"""
        return {
            'format': 'GFF3',
            'feature_types': self.feature_types,
            'strand_aware': True
        }


class dsRNADataLoader:
    """Handles loading and filtering of dsRNA data"""
    
    # Column name mappings: dsRNAscan TSV output -> internal names
    COLUMN_MAP_TSV = {
        'Chromosome': 'chr', 'Strand': 'strand',
        'i_start': 'i_start', 'i_end': 'i_end',
        'j_start': 'j_start', 'j_end': 'j_end',
        'eff_i_start': 'eff_i_start', 'eff_i_end': 'eff_i_end',
        'eff_j_start': 'eff_j_start', 'eff_j_end': 'eff_j_end',
        'stability_model_score': 'stability_model_score', 'probing_model_score': 'probing_model_score',
    }
    # Column name mappings: legacy parquet format -> internal names
    COLUMN_MAP_PARQUET = {
        'er_chr': 'chr', 'er_strand': 'strand',
        'er_i_start': 'i_start', 'er_i_end': 'i_end',
        'er_j_start': 'j_start', 'er_j_end': 'j_end',
    }

    def __init__(self, dsrna_file: Optional[str] = None):
        """Initialize with dsRNA data file (TSV from dsRNAscan or parquet)"""
        if dsrna_file:
            self.dsrna_file = dsrna_file
        else:
            # Try standard locations
            possible_paths = [
                "data/dsrna_predictions.parquet",
                os.environ.get("DSRNA_DATA_PATH", ""),
            ]
            possible_paths = [p for p in possible_paths if p]

            for path in possible_paths:
                if os.path.exists(path):
                    self.dsrna_file = path
                    break
            else:
                raise FileNotFoundError(
                    "dsRNA data file not found. Provide with --dsrna-file or\n"
                    "set DSRNA_DATA_PATH environment variable."
                )

        self.df = None
        self.available_subsets = {}
        self.file_format = None  # 'tsv' or 'parquet'
        
    def load_data(self, columns: Optional[List[str]] = None):
        """Load dsRNA data from TSV (dsRNAscan output) or parquet format.

        Automatically detects format and normalizes column names.
        """
        ext = Path(self.dsrna_file).suffix.lower()

        if ext in ['.tsv', '.txt', '.csv']:
            self._load_tsv()
        elif ext == '.parquet':
            self._load_parquet(columns)
        else:
            # Try TSV first (dsRNAscan output), fall back to parquet
            try:
                self._load_tsv()
            except Exception:
                self._load_parquet(columns)

        logger.info(f"Loaded {len(self.df):,} dsRNA regions from {self.file_format} file")
        self._define_subsets()

    def _load_tsv(self):
        """Load dsRNAscan TSV output and normalize columns"""
        sep = '\t'
        if self.dsrna_file.endswith('.csv'):
            sep = ','
        self.df = pd.read_csv(self.dsrna_file, sep=sep)
        self.file_format = 'tsv'

        # Rename columns to internal names
        rename = {k: v for k, v in self.COLUMN_MAP_TSV.items() if k in self.df.columns}
        self.df = self.df.rename(columns=rename)

        self.has_structure_probing = 'probing_model_score' in self.df.columns

    def _load_parquet(self, columns=None):
        """Load legacy parquet format and normalize columns"""
        if columns is None:
            columns = ['er_chr', 'er_i_start', 'er_i_end', 'er_j_start', 'er_j_end',
                       'er_strand']

        # Try loading with optional columns
        optional_cols = ['alu', 'i_phast17', 'j_phast17', 'i_phast100', 'j_phast100',
                        'pred_prob_editing_structure_only_alu100pct',
                        'pred_prob_editing_gtex_alu100pct',
                        'pred_3utr_normalized',
                        'pred_3utr_All_power_weighted_advantage']

        # Load with whatever columns exist
        all_cols = columns + optional_cols
        try:
            self.df = pd.read_parquet(self.dsrna_file, columns=all_cols)
        except Exception:
            # Fall back to just the required columns
            available = pd.read_parquet(self.dsrna_file, columns=None).columns.tolist()
            load_cols = [c for c in all_cols if c in available]
            self.df = pd.read_parquet(self.dsrna_file, columns=load_cols)

        self.file_format = 'parquet'

        # Rename columns to internal names
        rename = {k: v for k, v in self.COLUMN_MAP_PARQUET.items() if k in self.df.columns}
        self.df = self.df.rename(columns=rename)

        # Detect structure probing column
        if 'pred_3utr_normalized' in self.df.columns:
            self.has_structure_probing = True
            self.structure_probing_col = 'pred_3utr_normalized'
            self.structure_probing_threshold = 0.4574
        elif 'pred_3utr_All_power_weighted_advantage' in self.df.columns:
            self.has_structure_probing = True
            self.structure_probing_col = 'pred_3utr_All_power_weighted_advantage'
            self.structure_probing_threshold = 0.0315
        else:
            self.has_structure_probing = False
    
    def _define_subsets(self):
        """Define available dsRNA subsets based on columns present in the data"""
        self.available_subsets = {'all': None}

        # Conservation subsets (parquet format only)
        has_conservation = all(c in self.df.columns for c in
                              ['i_phast17', 'j_phast17', 'i_phast100', 'j_phast100'])
        if has_conservation:
            conserved_mask = (
                (self.df['i_phast17'] > 0.5) & (self.df['j_phast17'] > 0.5) &
                (self.df['i_phast100'] > 0.5) & (self.df['j_phast100'] > 0.5)
            )
            self.available_subsets['conserved'] = conserved_mask

        # Legacy ML masks (parquet format)
        ml_masks = []
        if 'pred_prob_editing_structure_only_alu100pct' in self.df.columns:
            ml_structure_mask = self.df['pred_prob_editing_structure_only_alu100pct'] > 0.2471
            self.available_subsets['ml_structure'] = ml_structure_mask
            ml_masks.append(ml_structure_mask)
        if 'pred_prob_editing_gtex_alu100pct' in self.df.columns:
            ml_gtex_mask = self.df['pred_prob_editing_gtex_alu100pct'] > 0.2513
            self.available_subsets['ml_gtex'] = ml_gtex_mask
            ml_masks.append(ml_gtex_mask)

        # dsRNAscan ML scores (TSV format)
        # Use pre-computed categorical columns if available, otherwise threshold
        if 'likely_edited' in self.df.columns:
            self.available_subsets['likely_edited'] = self.df['likely_edited']
            ml_masks.append(self.df['likely_edited'])
        elif 'stability_model_score' in self.df.columns:
            mask = self.df['stability_model_score'] > 0.2471
            self.available_subsets['likely_edited'] = mask
            ml_masks.append(mask)
        if 'likely_forms' in self.df.columns:
            self.available_subsets['likely_forms'] = self.df['likely_forms']
            ml_masks.append(self.df['likely_forms'])
        elif 'probing_model_score' in self.df.columns:
            mask = self.df['probing_model_score'] > 0.0315
            self.available_subsets['likely_forms'] = mask
            ml_masks.append(mask)

        # Structure probing (parquet format)
        if self.has_structure_probing and hasattr(self, 'structure_probing_col'):
            ml_probing_mask = self.df[self.structure_probing_col] > self.structure_probing_threshold
            self.available_subsets['ml_structure_probing'] = ml_probing_mask
            ml_masks.append(ml_probing_mask)

        # Combined high-confidence subset
        if ml_masks:
            ml_any_mask = ml_masks[0]
            for m in ml_masks[1:]:
                ml_any_mask = ml_any_mask | m
            self.available_subsets['any_high_conf'] = ml_any_mask
            if has_conservation:
                self.available_subsets['conserved_high_conf'] = conserved_mask & ml_any_mask

        # Alu-based subsets (parquet format)
        if 'alu' in self.df.columns:
            alu_mask = self.df['alu'] == 'Alu'
            nonalu_mask = self.df['alu'] == 'Non-Alu'
            self.available_subsets['alu'] = alu_mask
            self.available_subsets['nonalu'] = nonalu_mask
            if ml_masks:
                self.available_subsets['alu_high_conf'] = alu_mask & ml_any_mask
                self.available_subsets['nonalu_high_conf'] = nonalu_mask & ml_any_mask

        # Log available subsets
        for name, mask in self.available_subsets.items():
            count = len(self.df) if mask is None else mask.sum()
            logger.info(f"  Subset '{name}': {count:,} dsRNAs")
    
    def get_subset_bed(self, subset_name: str = 'all') -> Tuple[pybedtools.BedTool, int]:
        """Get BedTool for specified dsRNA subset (vectorized)"""
        if subset_name not in self.available_subsets:
            raise ValueError(f"Unknown subset: {subset_name}. Available: {list(self.available_subsets.keys())}")

        mask = self.available_subsets[subset_name]
        subset_df = self.df if mask is None else self.df[mask]

        logger.info(f"Creating BED entries for {len(subset_df):,} dsRNAs...")

        # Vectorized: build both arms as DataFrames and concatenate
        idx_labels = subset_df.index.astype(str)

        i_arms = pd.DataFrame({
            'chrom': subset_df['chr'].values,
            'start': subset_df['i_start'].astype(int).values,
            'end': subset_df['i_end'].astype(int).values,
            'name': 'dsRNA_' + idx_labels + '_i',
            'score': 0,
            'strand': subset_df['strand'].values,
        })
        j_arms = pd.DataFrame({
            'chrom': subset_df['chr'].values,
            'start': subset_df['j_start'].astype(int).values,
            'end': subset_df['j_end'].astype(int).values,
            'name': 'dsRNA_' + idx_labels + '_j',
            'score': 0,
            'strand': subset_df['strand'].values,
        })

        bed_df = pd.concat([i_arms, j_arms], ignore_index=True)

        temp_bed = get_temp_filename(f"dsrna_{subset_name}", ".bed")
        bed_df.to_csv(temp_bed, sep='\t', header=False, index=False)

        return pybedtools.BedTool(temp_bed), len(subset_df)


class StreamingOverlapAnalyzer:
    """Core overlap analysis module"""
    
    def __init__(self, query_bed: pybedtools.BedTool, dsrna_bed: pybedtools.BedTool, 
                 n_dsrnas: int, options: Dict):
        """Initialize analyzer with data"""
        self.query_bed = query_bed
        self.dsrna_bed = dsrna_bed
        self.n_dsrnas = n_dsrnas
        self.options = options
        self.results = {}
        self.control_beds = options.get('control_beds', [])
        
    def calculate_basic_overlaps(self):
        """Calculate basic overlap statistics"""
        # Count features
        n_query_features = len(self.query_bed)
        n_dsrna_intervals = len(self.dsrna_bed)
        
        # Calculate overlaps
        if self.options.get('strand_specific', False):
            # Strand-specific overlap
            query_with_dsrna = self.query_bed.intersect(self.dsrna_bed, u=True, s=True)
            dsrna_with_query = self.dsrna_bed.intersect(self.query_bed, u=True, s=True)
        else:
            # Strand-agnostic overlap
            query_with_dsrna = self.query_bed.intersect(self.dsrna_bed, u=True)
            dsrna_with_query = self.dsrna_bed.intersect(self.query_bed, u=True)
        
        n_query_overlapping = len(query_with_dsrna)
        n_dsrna_overlapping = len(dsrna_with_query)
        
        # Calculate percentages
        query_overlap_pct = 100 * n_query_overlapping / n_query_features if n_query_features > 0 else 0
        dsrna_overlap_pct = 100 * n_dsrna_overlapping / n_dsrna_intervals if n_dsrna_intervals > 0 else 0
        
        self.results['basic'] = {
            'n_query_features': n_query_features,
            'n_dsrna_regions': self.n_dsrnas,
            'n_dsrna_intervals': n_dsrna_intervals,
            'n_query_overlapping': n_query_overlapping,
            'n_dsrna_overlapping': n_dsrna_overlapping,
            'query_overlap_pct': query_overlap_pct,
            'dsrna_overlap_pct': dsrna_overlap_pct,
            'reciprocal_overlap_score': np.sqrt(query_overlap_pct * dsrna_overlap_pct / 100)
        }
        
        return self.results['basic']
    
    def calculate_enrichment(self, n_permutations: int = 100):
        """Calculate enrichment using permutation testing"""
        
        # Get observed overlap count
        observed_count = self.results['basic']['n_dsrna_overlapping']
        
        # Check if we have pre-generated control beds
        if self.control_beds:
            logger.info(f"Using {len(self.control_beds)} pre-generated control beds...")
            control_counts = self._use_pregenerated_controls()
        else:
            logger.info(f"Running {n_permutations} permutations for enrichment analysis...")
            control_counts = self._run_permutations(n_permutations)
        
        control_counts = np.array(control_counts)
        
        # Calculate statistics
        mean_control = np.mean(control_counts)
        std_control = np.std(control_counts)
        
        # Fold enrichment
        fold_enrichment = observed_count / mean_control if mean_control > 0 else np.inf
        
        # Z-score
        z_score = (observed_count - mean_control) / std_control if std_control > 0 else 0
        
        # Empirical p-values
        p_enrichment = (np.sum(control_counts >= observed_count) + 1) / (n_permutations + 1)
        p_depletion = (np.sum(control_counts <= observed_count) + 1) / (n_permutations + 1)
        
        # Two-tailed p-value
        p_value = 2 * min(p_enrichment, p_depletion)
        #p_value = min(p_value, 1.0)
        
        self.results['enrichment'] = {
            'observed_overlaps': observed_count,
            'expected_overlaps': mean_control,
            'control_mean': mean_control,
            'control_std': std_control,
            'control_min': np.min(control_counts),
            'control_max': np.max(control_counts),
            'fold_enrichment': fold_enrichment,
            'z_score': z_score,
            'p_value': p_value,
            'p_enrichment': p_enrichment,
            'p_depletion': p_depletion,
            'n_permutations': n_permutations
        }
        
        return self.results['enrichment']
    
    def _get_chromosome_sizes(self) -> Dict[str, int]:
        """Get chromosome sizes from the dsRNA data"""
        chrom_sizes = {}
        
        # Get unique chromosomes and their max coordinates
        for interval in self.dsrna_bed:
            chrom = interval.chrom
            end = interval.end
            
            if chrom not in chrom_sizes:
                chrom_sizes[chrom] = end
            else:
                chrom_sizes[chrom] = max(chrom_sizes[chrom], end)
        
        # Add some buffer to ensure we don't place intervals at the very end
        for chrom in chrom_sizes:
            chrom_sizes[chrom] = int(chrom_sizes[chrom] * 1.1)
        
        return chrom_sizes
    
    def _run_permutations(self, n_permutations: int) -> List[int]:
        """Run permutation testing with streaming approach to minimize temp files"""
        # Get chromosome sizes for shuffling
        chrom_sizes = self._get_chromosome_sizes()
        
        # Run permutations
        control_counts = []
        
        # Create a single reusable temp file
        temp_bed_file = get_temp_filename("shuffled", ".bed")
        
        for i in tqdm(range(n_permutations), desc="Running permutations"):
            # Stream shuffled regions directly to file
            with open(temp_bed_file, 'w') as f:
                for interval in self.dsrna_bed:
                    chrom = interval.chrom
                    length = interval.end - interval.start
                    strand = interval.strand if len(interval.fields) >= 6 else '.'
                    
                    if chrom in chrom_sizes:
                        max_pos = chrom_sizes[chrom] - length
                        if max_pos > 0:
                            new_start = np.random.randint(0, max_pos)
                            new_end = new_start + length
                            f.write(f"{chrom}\t{new_start}\t{new_end}\tshuffled\t0\t{strand}\n")
            
            # Create BedTool from the file
            shuffled_dsrna = pybedtools.BedTool(temp_bed_file)
            
            # Count overlaps with shuffled regions
            if self.options.get('strand_specific', False):
                control_overlap = len(shuffled_dsrna.intersect(self.query_bed, u=True, s=True))
            else:
                control_overlap = len(shuffled_dsrna.intersect(self.query_bed, u=True))
            
            control_counts.append(control_overlap)
            
        return control_counts
    
    def _use_pregenerated_controls(self) -> List[int]:
        """Use pre-generated control BED files for faster analysis"""
        control_counts = []
        
        for control_file in tqdm(self.control_beds, desc="Processing control files"):
            try:
                control_bed = pybedtools.BedTool(control_file)
                
                # Count overlaps with control regions
                if self.options.get('strand_specific', False):
                    control_overlap = len(control_bed.intersect(self.query_bed, u=True, s=True))
                else:
                    control_overlap = len(control_bed.intersect(self.query_bed, u=True))
                
                control_counts.append(control_overlap)
            except Exception as e:
                logger.warning(f"Error processing control file {control_file}: {e}")
                continue
        
        if len(control_counts) < 10:
            logger.warning(f"Only {len(control_counts)} valid control files found. Results may be unreliable.")
        
        return control_counts
    


def detect_file_format(filepath: str) -> str:
    """Detect file format from extension and content"""
    ext = Path(filepath).suffix.lower()
    
    if ext in ['.bed', '.narrowpeak', '.broadpeak']:
        return 'BED'
    elif ext in ['.gff', '.gff3']:
        return 'GFF3'
    elif ext in ['.gtf']:
        return 'GTF'
    else:
        # Try to detect from content
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    try:
                        int(parts[1])
                        int(parts[2])
                        return 'BED'
                    except:
                        if len(parts) >= 9:
                            return 'GFF3'
                break
    
    return 'UNKNOWN'


def main():
    """Main analysis function"""
    # Setup temp directory at the start
    setup_temp_dir()
    
    parser = argparse.ArgumentParser(
        description="Analyze overlap between genomic features and dsRNA regions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Available dsRNA subsets:
  all                  - All dsRNA regions (5.1M, no filtering)
  conserved            - PhastCons > 0.5 in both arms (10.7K)
  conserved_high_conf  - Conserved AND ML high confidence (1.6K) 
  ml_structure         - ML Structure-only model high confidence (2.1M)
  ml_gtex              - ML GTEx model high confidence (1.7M)
  ml_structure_probing - Structure probing validated (2.1M)
  any_high_conf        - Any ML model high confidence (2.4M)
  alu                  - Alu-derived dsRNAs (2.6M)
  nonalu               - Non-Alu dsRNAs (2.5M)
  alu_high_conf        - Alu AND ML high confidence (2.1M)
  nonalu_high_conf     - Non-Alu AND ML high confidence (310K)

Recommended subsets by dataset size:
  Small (<100 features)      : conserved_high_conf, conserved
  Medium (100-1000 features) : conserved, ml_structure, nonalu_high_conf
  Large (1000-10000 features): any_high_conf, ml_gtex, alu/nonalu separate
  Very Large (>10000)        : all, any_high_conf

Examples:
  # Basic BED file analysis (enhanced output and plots are now default)
  python dsrna_overlap_analyzer.py my_peaks.bed
  
  # GFF3 with specific features
  python dsrna_overlap_analyzer.py annotations.gff3 --feature-type gene,exon
  
  # High-confidence conserved dsRNAs with enrichment analysis
  python dsrna_overlap_analyzer.py chip_seq.bed --subset conserved_high_conf --permutations 100
  
  # Strand-specific analysis with simple output (no colors/plots)
  python dsrna_overlap_analyzer.py rbp_clips.bed --strand-specific --no-enhanced-output --no-plot
  
  # Quick test with subset recommendation
  python dsrna_overlap_analyzer.py my_data.bed --recommend-subset
        """
    )
    
    # Required arguments (unless using special modes)
    parser.add_argument('input_file', nargs='?', help='Input file (BED/GFF3/GTF format)')
    
    # Optional arguments
    parser.add_argument('--format', choices=['BED', 'GFF3', 'GTF', 'auto'], 
                       default='auto', help='Input file format (default: auto-detect)')
    
    parser.add_argument('--subset', default='conserved_high_conf',
                       help='dsRNA subset to analyze (see available subsets above, default: conserved_high_conf)')

    parser.add_argument('--feature-type', nargs='+',
                       help='For GFF3/GTF: specific feature types to analyze')
    
    parser.add_argument('--strand-specific', action='store_true',
                       help='Perform strand-specific overlap analysis')
    
    parser.add_argument('--min-overlap', type=float, default=1,
                       help='Minimum overlap in base pairs (default: 1)')
    
    parser.add_argument('--permutations', type=int, default=100,
                       help='Number of permutations for enrichment testing (default: 100)')

    parser.add_argument('--control-dir', type=str,
                       help='Directory containing pre-generated control BED files (faster than permutations)')
    
    parser.add_argument('--output-prefix', default='dsrna_overlap',
                       help='Prefix for output files (default: dsrna_overlap)')
    
    parser.add_argument('--output-format', choices=['json', 'csv', 'both'], default='both',
                       help='Output format (default: both - CSV for viewing, JSON for metadata)')
    
    parser.add_argument('--dsrna-file', 
                       help='Custom dsRNA data file (default: use standard file)')
    
    parser.add_argument('--list-subsets', action='store_true',
                       help='List available dsRNA subsets and exit')
    
    parser.add_argument('--generate-controls', type=int,
                       help='Generate N control BED files for the specified subset and exit')
    
    parser.add_argument('--recommend-subset', action='store_true',
                       help='Recommend optimal subset based on input file size')
    
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose output')
    
    parser.add_argument('--keep-temp', action='store_true',
                       help='Keep temporary files (for debugging)')
    
    parser.add_argument('--enhanced-output', action='store_true', default=True,
                       help='Use enhanced colorful output with biological interpretation (default: enabled)')
    
    parser.add_argument('--no-enhanced-output', dest='enhanced_output', action='store_false',
                       help='Disable enhanced output and use simple format')
    
    parser.add_argument('--plot', action='store_true', default=True,
                       help='Generate enrichment visualization plot (default: enabled)')
    
    parser.add_argument('--no-plot', dest='plot', action='store_false',
                       help='Disable plot generation')
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # If keeping temp files, disable cleanup
    if args.keep_temp:
        atexit.unregister(cleanup_temp_dir)
        logger.info(f"Temp files will be kept in: {TEMP_DIR}")
    
    # Initialize dsRNA data loader
    logger.info("Loading dsRNA data...")
    dsrna_loader = dsRNADataLoader(args.dsrna_file)
    dsrna_loader.load_data()
    
    # List subsets if requested
    if args.list_subsets:
        print("\nAvailable dsRNA subsets:")
        print("-" * 60)
        print(f"{'Subset':<25} {'dsRNAs':>12} {'Description'}")
        print("-" * 60)
        
        subset_descriptions = {
            'all': 'All dsRNA regions (no filtering)',
            'conserved': 'PhastCons > 0.5 in both arms',
            'ml_structure': 'ML Structure-only model high confidence',
            'ml_gtex': 'ML GTEx model high confidence',
            'ml_structure_probing': 'Structure probing validated',
            'any_high_conf': 'Any ML model high confidence',
            'conserved_high_conf': 'Conserved AND ML high confidence',
            'alu': 'Alu-derived dsRNAs',
            'nonalu': 'Non-Alu dsRNAs',
            'alu_high_conf': 'Alu AND ML high confidence',
            'nonalu_high_conf': 'Non-Alu AND ML high confidence'
        }
        
        for subset, mask in dsrna_loader.available_subsets.items():
            if mask is None:
                count = len(dsrna_loader.df)
            else:
                count = mask.sum()
            desc = subset_descriptions.get(subset, '')
            print(f"{subset:<25} {count:>12,} {desc}")
        
        print("\nRecommended subsets by use case:")
        print("  - Quick test: 'conserved' or 'conserved_high_conf'")
        print("  - Standard analysis: 'any_high_conf' or 'ml_structure'")
        print("  - Alu-specific: 'alu_high_conf' or 'nonalu_high_conf'")
        print("  - Comprehensive: 'all' (warning: very large)")
        sys.exit(0)
    
    # Generate controls if requested
    if args.generate_controls:
        if not args.subset:
            print("Error: --subset must be specified with --generate-controls")
            sys.exit(1)
        
        print(f"Generating {args.generate_controls} control files for subset '{args.subset}'...")
        dsrna_bed, n_dsrnas = dsrna_loader.get_subset_bed(args.subset)
        
        # Create output directory
        control_dir = Path(f"controls_{args.subset}")
        control_dir.mkdir(exist_ok=True)
        
        # Get chromosome sizes
        chrom_sizes = {}
        for interval in dsrna_bed:
            chrom = interval.chrom
            end = interval.end
            if chrom not in chrom_sizes:
                chrom_sizes[chrom] = end
            else:
                chrom_sizes[chrom] = max(chrom_sizes[chrom], end)
        
        # Add buffer
        for chrom in chrom_sizes:
            chrom_sizes[chrom] = int(chrom_sizes[chrom] * 1.1)
        
        # Generate control files
        from tqdm import tqdm
        for i in tqdm(range(args.generate_controls), desc="Generating controls"):
            # Create shuffled regions
            shuffled_entries = []
            for interval in dsrna_bed:
                chrom = interval.chrom
                length = interval.end - interval.start
                strand = interval.strand if len(interval.fields) >= 6 else '.'
                name = f"control_{i}_{interval.name}"
                
                if chrom in chrom_sizes:
                    max_pos = chrom_sizes[chrom] - length
                    if max_pos > 0:
                        new_start = np.random.randint(0, max_pos)
                        new_end = new_start + length
                        shuffled_entries.append([chrom, new_start, new_end, name, 0, strand])
            
            # Save control file
            control_file = control_dir / f"control_{i:04d}.bed"
            pd.DataFrame(shuffled_entries).to_csv(control_file, sep='\t', header=False, index=False)
            
            # Compress
            os.system(f"gzip {control_file}")
        
        print(f"\nGenerated {args.generate_controls} control files in {control_dir}/")
        print(f"Total size: {sum(f.stat().st_size for f in control_dir.glob('*.gz')) / 1024 / 1024:.1f} MB")
        sys.exit(0)
    
    # Check if we need input file
    if not args.input_file:
        logger.error("Input file is required for analysis")
        parser.print_help()
        sys.exit(1)
    
    # Check input file
    if not os.path.exists(args.input_file):
        logger.error(f"Input file not found: {args.input_file}")
        sys.exit(1)
    
    # Detect format
    if args.format == 'auto':
        file_format = detect_file_format(args.input_file)
        logger.info(f"Detected file format: {file_format}")
    else:
        file_format = args.format
    
    # Parse input file
    logger.info(f"Parsing {file_format} file: {args.input_file}")
    
    if file_format == 'BED':
        parser = BEDParser()
    elif file_format == 'GFF3':
        parser = GFF3Parser(feature_types=args.feature_type)
    else:
        logger.error(f"Unsupported format: {file_format}")
        sys.exit(1)
    
    query_bed = parser.parse(args.input_file)
    metadata = parser.get_metadata()
    
    # Get dsRNA subset
    logger.info(f"Loading dsRNA subset: {args.subset}")
    dsrna_bed, n_dsrnas = dsrna_loader.get_subset_bed(args.subset)
    
    # Load control beds if provided
    control_beds = []
    if args.control_dir:
        control_dir = Path(args.control_dir)
        if control_dir.exists():
            # Look for control files specific to this subset
            subset_controls = list(control_dir.glob(f"*{args.subset}*.bed*"))
            if not subset_controls:
                # Fall back to general controls
                subset_controls = list(control_dir.glob("randomized_*.bed*"))
            
            control_beds = [str(f) for f in subset_controls[:100]]  # Use up to 100 controls
            logger.info(f"Found {len(control_beds)} control BED files")
        else:
            logger.warning(f"Control directory not found: {args.control_dir}")
    
    # Analyze overlaps
    logger.info("Calculating overlaps...")
    analyzer = StreamingOverlapAnalyzer(
        query_bed, 
        dsrna_bed, 
        n_dsrnas,
        {
            'strand_specific': args.strand_specific,
            'min_overlap': args.min_overlap,
            'control_beds': control_beds
        }
    )
    
    # Basic overlap analysis
    basic_results = analyzer.calculate_basic_overlaps()
    basic_results['format'] = file_format  # Add format to results
    
    # Enrichment analysis if requested
    enrichment_results = None
    if args.permutations > 0 or control_beds:
        n_perms = args.permutations if not control_beds else len(control_beds)
        enrichment_results = analyzer.calculate_enrichment(n_perms)
    
    # Generate timestamp for output files
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    # Print results using enhanced output if available and requested
    if args.enhanced_output and HAS_ENHANCED_OUTPUT:
        print_enhanced_results(basic_results, enrichment_results, args)
        
        # Generate plot if requested
        if args.plot and enrichment_results:
            plot_file = f"{args.output_prefix}_enrichment_plot_{timestamp}.png"
            create_enrichment_plot(enrichment_results, plot_file)
        
        # Print subset recommendation if requested
        if args.recommend_subset:
            print_subset_recommendation(basic_results['n_query_features'])
    else:
        # Fall back to original output format
        print("\n" + "="*60)
        print("dsRNA OVERLAP ANALYSIS RESULTS")
        print("="*60)
        print(f"Input file: {args.input_file}")
        print(f"Format: {file_format}")
        print(f"dsRNA subset: {args.subset}")
        print(f"Strand-specific: {args.strand_specific}")
        print("-"*60)
        print(f"Query features: {basic_results['n_query_features']:,}")
        print(f"dsRNA regions: {basic_results['n_dsrna_regions']:,}")
        print(f"dsRNA intervals (i+j arms): {basic_results['n_dsrna_intervals']:,}")
        print("-"*60)
        print(f"Query features overlapping dsRNA: {basic_results['n_query_overlapping']:,} ({basic_results['query_overlap_pct']:.1f}%)")
        print(f"dsRNA intervals overlapping query: {basic_results['n_dsrna_overlapping']:,} ({basic_results['dsrna_overlap_pct']:.1f}%)")
        print(f"Reciprocal overlap score: {basic_results['reciprocal_overlap_score']:.3f}")
        
        if enrichment_results:
            print("-"*60)
            print("ENRICHMENT ANALYSIS (Permutation Test)")
            print("-"*60)
            print(f"Observed overlaps: {enrichment_results['observed_overlaps']:,}")
            print(f"Expected overlaps: {enrichment_results['expected_overlaps']:.1f} ± {enrichment_results['control_std']:.1f}")
            print(f"Control range: {enrichment_results['control_min']} - {enrichment_results['control_max']}")
            print(f"Fold enrichment: {enrichment_results['fold_enrichment']:.2f}x")
            print(f"Z-score: {enrichment_results['z_score']:.2f}")
            print(f"P-value (two-tailed): {enrichment_results['p_value']:.4f}")
            
            # Interpretation
            if enrichment_results['p_value'] < 0.05:
                if enrichment_results['fold_enrichment'] > 1:
                    print(f"\n✓ SIGNIFICANT ENRICHMENT: {enrichment_results['fold_enrichment']:.2f}x more overlap than expected by chance")
                else:
                    print(f"\n✓ SIGNIFICANT DEPLETION: {enrichment_results['fold_enrichment']:.2f}x less overlap than expected by chance")
            else:
                print("\n✗ No significant enrichment or depletion detected")
    
    # Save detailed results
    # (timestamp already generated above)
    
    # Convert numpy types to Python types for JSON serialization
    def convert_to_serializable(obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert_to_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_to_serializable(v) for v in obj]
        return obj
    
    results = {
        'metadata': {
            'input_file': args.input_file,
            'format': file_format,
            'subset': args.subset,
            'timestamp': timestamp,
            'options': vars(args)
        },
        'results': convert_to_serializable(analyzer.results)
    }
    
    # Save results based on output format
    if args.output_format in ['json', 'both']:
        json_file = f"{args.output_prefix}_results_{timestamp}.json"
        with open(json_file, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\n📊 JSON results saved to: {json_file} (complete metadata & detailed results)")
    
    if args.output_format in ['csv', 'both']:
        csv_file = f"{args.output_prefix}_results_{timestamp}.csv"
        
        # Create CSV-friendly summary
        csv_data = {
            'input_file': args.input_file,
            'feature_name': Path(args.input_file).stem,
            'dsrna_subset': args.subset,
            'total_features': basic_results['n_query_features'],
            'total_dsrna_regions': basic_results['n_dsrna_regions'],
            'features_overlapping_dsrna': basic_results['n_query_overlapping'],
            'dsrna_overlapping_features': basic_results['n_dsrna_overlapping'],
            'feature_overlap_percent': round(basic_results['query_overlap_pct'], 2),
            'dsrna_overlap_percent': round(basic_results['dsrna_overlap_pct'], 2),
            'reciprocal_overlap_score': round(basic_results['reciprocal_overlap_score'], 3)
        }
        
        # Add enrichment data if available
        if 'enrichment' in analyzer.results:
            enrich = analyzer.results['enrichment']
            csv_data.update({
                'observed_overlaps': enrich['observed_overlaps'],
                'expected_overlaps': round(enrich['expected_overlaps'], 1),
                'fold_enrichment': round(enrich['fold_enrichment'], 2),
                'z_score': round(enrich['z_score'], 2),
                'p_value': round(enrich['p_value'], 4),
                'n_permutations': enrich['n_permutations']
            })
        
        # Write CSV
        with open(csv_file, 'w') as f:
            # Write header
            f.write(','.join(csv_data.keys()) + '\n')
            # Write values
            f.write(','.join(str(v) for v in csv_data.values()) + '\n')
        
        print(f"📁 CSV results saved to: {csv_file} (summary table for Excel/R)")
    
    if not args.keep_temp:
        print(f"Temporary files cleaned up automatically")
    else:
        print(f"Temporary files kept in: {TEMP_DIR}")
    
    # pybedtools cleanup (removes its own temp files)
    pybedtools.cleanup()


if __name__ == "__main__":
    main()