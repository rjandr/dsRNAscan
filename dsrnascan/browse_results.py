#!/usr/bin/env python3
"""
dsRNAscan Results Browser - Simple Forna Viewer

Quick visualization of dsRNA structures from dsRNAscan output.
Supports optional RNA editing site annotation from BED or GFF3 files.
"""

import os
import sys
import json
import pandas as pd
import numpy as np
import socketserver
import webbrowser
import argparse
from pathlib import Path
from http.server import SimpleHTTPRequestHandler
from urllib.parse import urlparse, parse_qs
from collections import defaultdict
import bisect

class CORSHTTPRequestHandler(SimpleHTTPRequestHandler):
    """HTTP request handler with CORS support"""
    
    def end_headers(self):
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type')
        super().end_headers()
    
    def do_OPTIONS(self):
        self.send_response(200)
        self.end_headers()
    
    def handle_error(self, request, client_address):
        """Override to suppress broken pipe errors"""
        import errno
        exc_type, exc_value = sys.exc_info()[:2]
        if exc_type == BrokenPipeError or (exc_type == OSError and exc_value.errno == errno.EPIPE):
            # Ignore broken pipe errors
            pass
        else:
            super().handle_error(request, client_address)

class OptimizedEditingSites:
    """Optimized storage and lookup for editing sites using binary search"""
    
    def __init__(self):
        # Store as sorted arrays per chromosome/strand for binary search
        # Using numpy arrays for memory efficiency
        self.sites = defaultdict(lambda: defaultdict(lambda: {'positions': [], 'frequencies': []}))
        self.finalized = False
    
    def add_site(self, chrom, strand, position, frequency):
        """Add a site while maintaining sorted order"""
        if self.finalized:
            raise RuntimeError("Cannot add sites after finalization")
        
        sites = self.sites[chrom][strand]
        # Use bisect to maintain sorted order during insertion
        idx = bisect.bisect_left(sites['positions'], position)
        sites['positions'].insert(idx, position)
        sites['frequencies'].insert(idx, frequency)
    
    def finalize(self):
        """Convert lists to numpy arrays for memory efficiency"""
        print("Optimizing editing site storage...")
        for chrom in self.sites:
            for strand in self.sites[chrom]:
                data = self.sites[chrom][strand]
                if data['positions']:
                    # Use uint32 for positions (saves memory for large files)
                    data['positions'] = np.array(data['positions'], dtype=np.uint32)
                    # Use float16 for frequencies (sufficient precision, saves memory)
                    data['frequencies'] = np.array(data['frequencies'], dtype=np.float16)
        self.finalized = True
    
    def find_in_range(self, chrom, strand, start, end):
        """Binary search for sites in range - O(log n) instead of O(n)"""
        if chrom not in self.sites or strand not in self.sites[chrom]:
            return []
        
        sites = self.sites[chrom][strand]
        positions = sites['positions']
        
        if len(positions) == 0:
            return []
        
        # Binary search for range
        if self.finalized:
            # NumPy searchsorted for arrays
            start_idx = np.searchsorted(positions, start, side='left')
            end_idx = np.searchsorted(positions, end, side='right')
        else:
            # bisect for lists
            start_idx = bisect.bisect_left(positions, start)
            end_idx = bisect.bisect_right(positions, end)
        
        return [(int(positions[i]), float(sites['frequencies'][i])) 
                for i in range(start_idx, end_idx)]
    
    def get_summary(self):
        """Get summary statistics"""
        total_sites = 0
        chroms = set()
        for chrom in self.sites:
            chroms.add(chrom)
            for strand in self.sites[chrom]:
                total_sites += len(self.sites[chrom][strand]['positions'])
        return {'total_sites': total_sites, 'chromosomes': len(chroms)}

def parse_editing_file_optimized(editing_file, needed_chromosomes=None):
    """
    Parse BED or GFF3 file with editing sites using optimized storage.
    
    Args:
        editing_file: Path to BED or GFF3 file
        needed_chromosomes: Set of chromosomes to load (None = load all)
    
    Returns:
        OptimizedEditingSites object
    """
    editing_sites = OptimizedEditingSites()
    
    # First, count total lines for progress reporting
    print(f"Counting entries in {editing_file}...")
    total_lines = 0
    with open(editing_file, 'r') as f:
        for line in f:
            if not line.startswith('#') and line.strip():
                total_lines += 1
    
    if total_lines > 10000:
        print(f"Large file detected: {total_lines:,} entries. Processing...")
    
    try:
        # Detect file format
        is_gff3 = False
        with open(editing_file, 'r') as f:
            for line in f:
                if line.startswith('##gff-version'):
                    is_gff3 = True
                    break
                elif not line.startswith('#') and line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) == 9 and parts[2] not in ['', '.']:
                        is_gff3 = True
                    break
        
        # Parse based on format
        processed = 0
        skipped = 0
        
        with open(editing_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                chrom = parts[0]
                
                # Skip chromosomes we don't need (lazy loading)
                if needed_chromosomes and chrom not in needed_chromosomes:
                    skipped += 1
                    processed += 1
                    continue
                
                if is_gff3:
                    # GFF3 format parsing
                    if len(parts) < 9:
                        continue
                    
                    feature_type = parts[2]
                    
                    # Filter for editing-related features
                    if not any(term in feature_type.lower() for term in ['edit', 'modification', 'variant', 'snp', 'snv']):
                        attributes = parts[8]
                        if not any(term in attributes.lower() for term in ['edit', 'adar', 'apobec', 'a-to-i', 'c-to-u']):
                            skipped += 1
                            processed += 1
                            continue
                    
                    position = int(parts[3])  # GFF3 is 1-based
                    strand = parts[6]
                    
                    # Get frequency from score or attributes
                    frequency = 1.0
                    if parts[5] != '.':
                        try:
                            score = float(parts[5])
                            if score <= 1.0:
                                frequency = score
                            elif score <= 100:
                                frequency = score / 100.0
                            else:
                                frequency = score / 1000.0
                        except:
                            pass
                    
                    # Parse attributes for frequency
                    attributes = parts[8]
                    for attr in attributes.split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            key = key.strip().lower()
                            if key in ['frequency', 'freq', 'confidence', 'score', 'editing_level']:
                                try:
                                    freq_val = float(value.strip('%'))
                                    if freq_val > 1.0:
                                        frequency = freq_val / 100.0
                                    else:
                                        frequency = freq_val
                                    break
                                except:
                                    pass
                
                else:
                    # BED format parsing
                    if len(parts) < 6:
                        continue
                    
                    position = int(parts[1]) + 1  # BED is 0-based, convert to 1-based
                    strand = parts[5]
                    
                    frequency = 1.0
                    if len(parts) > 4:
                        try:
                            score = float(parts[4])
                            if score <= 1.0:
                                frequency = score
                            else:
                                frequency = score / 1000.0
                        except:
                            pass
                
                # Add to optimized structure
                editing_sites.add_site(chrom, strand, position, frequency)
                
                processed += 1
                # Show progress for large files
                if total_lines > 10000 and processed % 10000 == 0:
                    pct = (processed / total_lines) * 100
                    print(f"  Processed {processed:,} / {total_lines:,} ({pct:.1f}%)", end='\r')
        
        # Clear progress line
        if total_lines > 10000:
            print(f"\n  Loaded {processed - skipped:,} sites, skipped {skipped:,}")
        
        # Finalize for memory efficiency
        editing_sites.finalize()
        
        # Report summary
        summary = editing_sites.get_summary()
        format_name = "GFF3" if is_gff3 else "BED"
        print(f"Loaded {summary['total_sites']:,} editing sites from {summary['chromosomes']} chromosomes ({format_name} format)")
        
        return editing_sites
    
    except Exception as e:
        print(f"Warning: Could not parse editing file: {e}")
        return OptimizedEditingSites()

def parse_editing_file(editing_file):
    """
    Parse BED or GFF3 file with editing sites.
    
    BED format: chr start end name score strand [frequency]
    GFF3 format: chr source type start end score strand phase attributes
    
    Returns dict: {chr: {strand: [(position, frequency)]}}
    """
    editing_sites = defaultdict(lambda: defaultdict(list))
    
    try:
        # Detect file format
        is_gff3 = False
        with open(editing_file, 'r') as f:
            for line in f:
                if line.startswith('##gff-version'):
                    is_gff3 = True
                    break
                elif not line.startswith('#') and line.strip():
                    # Check if it looks like GFF3 (9 tab-separated fields)
                    parts = line.strip().split('\t')
                    if len(parts) == 9 and parts[2] not in ['', '.']:
                        is_gff3 = True
                    break
        
        # Parse based on format
        with open(editing_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                
                if is_gff3:
                    # GFF3 format parsing
                    if len(parts) < 9:
                        continue
                    
                    chrom = parts[0]
                    feature_type = parts[2]
                    
                    # Filter for editing-related features
                    if not any(term in feature_type.lower() for term in ['edit', 'modification', 'variant', 'snp', 'snv']):
                        # Also check attributes for editing-related terms
                        attributes = parts[8]
                        if not any(term in attributes.lower() for term in ['edit', 'adar', 'apobec', 'a-to-i', 'c-to-u']):
                            continue
                    
                    # GFF3 is 1-based
                    position = int(parts[3])
                    strand = parts[6]
                    
                    # Try to get frequency from score or attributes
                    frequency = 1.0  # Default frequency
                    
                    # Check score field
                    if parts[5] != '.':
                        try:
                            score = float(parts[5])
                            if score <= 1.0:
                                frequency = score
                            elif score <= 100:
                                frequency = score / 100.0
                            else:
                                frequency = score / 1000.0
                        except:
                            pass
                    
                    # Parse attributes for frequency/confidence
                    attributes = parts[8]
                    for attr in attributes.split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            key = key.strip().lower()
                            if key in ['frequency', 'freq', 'confidence', 'score', 'editing_level']:
                                try:
                                    freq_val = float(value.strip('%'))
                                    if freq_val > 1.0:
                                        frequency = freq_val / 100.0
                                    else:
                                        frequency = freq_val
                                    break
                                except:
                                    pass
                    
                else:
                    # BED format parsing
                    if len(parts) < 6:
                        continue
                    
                    chrom = parts[0]
                    # BED is 0-based, convert to 1-based for genomic coords
                    position = int(parts[1]) + 1  
                    strand = parts[5]
                    
                    # Try to get frequency from score or additional column
                    frequency = 1.0  # Default frequency
                    if len(parts) > 4:
                        try:
                            score = float(parts[4])
                            if score <= 1.0:
                                frequency = score
                            else:
                                frequency = score / 1000.0
                        except:
                            pass
                
                # Store as 1-based genomic position
                editing_sites[chrom][strand].append((position, frequency))
        
        # Sort positions for each chromosome/strand
        for chrom in editing_sites:
            for strand in editing_sites[chrom]:
                editing_sites[chrom][strand].sort(key=lambda x: x[0])
        
        # Report what was loaded
        total_sites = sum(len(sites) for chrom in editing_sites for sites in editing_sites[chrom].values())
        format_name = "GFF3" if is_gff3 else "BED"
        print(f"Loaded {total_sites} editing sites from {format_name} file")
        
        return editing_sites
    
    except Exception as e:
        print(f"Warning: Could not parse editing file: {e}")
        return {}

def map_editing_to_dsrna_optimized(dsrna_row, editing_sites):
    """
    Map editing sites to dsRNA structure positions using optimized binary search.
    
    For forward strand: positions map directly
    For reverse strand: need to reverse the mapping
    
    Returns list of [position_in_structure, frequency] pairs
    """
    if not editing_sites:
        return []
    
    chrom = str(dsrna_row['Chromosome'])
    strand = str(dsrna_row['Strand'])
    
    # Get the effective coordinates (trimmed regions)
    if 'eff_i_start' in dsrna_row and pd.notna(dsrna_row['eff_i_start']):
        i_start = int(dsrna_row['eff_i_start'])
        i_end = int(dsrna_row['eff_i_end'])
        j_start = int(dsrna_row['eff_j_start'])
        j_end = int(dsrna_row['eff_j_end'])
    else:
        # Fall back to original coordinates
        i_start = int(dsrna_row['i_start'])
        i_end = int(dsrna_row['i_end'])
        j_start = int(dsrna_row['j_start'])
        j_end = int(dsrna_row['j_end'])
    
    # Get sequences to determine lengths
    i_seq = str(dsrna_row['i_seq']) if pd.notna(dsrna_row['i_seq']) else ''
    j_seq = str(dsrna_row['j_seq']) if pd.notna(dsrna_row['j_seq']) else ''
    i_length = len(i_seq)
    j_length = len(j_seq)
    
    structure_edits = []
    
    # Use optimized binary search to find sites in i-arm range
    if strand == '+':
        # Forward strand
        i_sites = editing_sites.find_in_range(chrom, strand, i_start, i_end)
        for edit_pos, frequency in i_sites:
            structure_pos = edit_pos - i_start
            if 0 <= structure_pos < i_length:
                structure_edits.append([int(structure_pos), float(frequency)])
        
        # Find sites in j-arm range
        j_sites = editing_sites.find_in_range(chrom, strand, j_start, j_end)
        for edit_pos, frequency in j_sites:
            structure_pos = i_length + (edit_pos - j_start)
            if i_length <= structure_pos < (i_length + j_length):
                structure_edits.append([int(structure_pos), float(frequency)])
    
    else:  # strand == '-'
        # Reverse strand: coordinates go from high to low (3' to 5')
        # i_start is at the 3' end, i_end is at the 5' end (i_end < i_start)
        i_sites = editing_sites.find_in_range(chrom, strand, i_end, i_start)
        for edit_pos, frequency in i_sites:
            pos_from_3prime = edit_pos - i_end
            structure_pos = i_length - 1 - pos_from_3prime
            if 0 <= structure_pos < i_length:
                structure_edits.append([int(structure_pos), float(frequency)])
        
        # j-arm (j_start < j_end for minus)
        j_sites = editing_sites.find_in_range(chrom, strand, j_start, j_end)
        for edit_pos, frequency in j_sites:
            pos_from_3prime = edit_pos - j_start
            structure_pos = i_length + (j_length - 1 - pos_from_3prime)
            if i_length <= structure_pos < (i_length + j_length):
                structure_edits.append([int(structure_pos), float(frequency)])
    
    return structure_edits

def map_editing_to_dsrna(dsrna_row, editing_sites):
    """
    Map editing sites to dsRNA structure positions (legacy version for compatibility).
    
    For forward strand: positions map directly
    For reverse strand: need to reverse the mapping
    
    Returns list of [position_in_structure, frequency] pairs
    """
    # Check if we have the optimized structure
    if isinstance(editing_sites, OptimizedEditingSites):
        return map_editing_to_dsrna_optimized(dsrna_row, editing_sites)
    
    # Legacy code for backward compatibility
    if not editing_sites:
        return []
    
    chrom = str(dsrna_row['Chromosome'])
    strand = str(dsrna_row['Strand'])
    
    # Get editing sites for this chromosome/strand
    if chrom not in editing_sites or strand not in editing_sites[chrom]:
        return []
    
    chrom_edits = editing_sites[chrom][strand]
    
    # Get the effective coordinates (trimmed regions)
    if 'eff_i_start' in dsrna_row and pd.notna(dsrna_row['eff_i_start']):
        i_start = int(dsrna_row['eff_i_start'])
        i_end = int(dsrna_row['eff_i_end'])
        j_start = int(dsrna_row['eff_j_start'])
        j_end = int(dsrna_row['eff_j_end'])
    else:
        # Fall back to original coordinates
        i_start = int(dsrna_row['i_start'])
        i_end = int(dsrna_row['i_end'])
        j_start = int(dsrna_row['j_start'])
        j_end = int(dsrna_row['j_end'])
    
    structure_edits = []
    
    # Get sequences to determine lengths
    i_seq = str(dsrna_row['i_seq']) if pd.notna(dsrna_row['i_seq']) else ''
    j_seq = str(dsrna_row['j_seq']) if pd.notna(dsrna_row['j_seq']) else ''
    i_length = len(i_seq)
    j_length = len(j_seq)
    
    # Map editing sites to structure positions
    for edit_pos, frequency in chrom_edits:
        structure_pos = None
        
        if strand == '+':
            # Forward strand: direct mapping
            # Check if edit is in i-arm
            if i_start <= edit_pos <= i_end:
                # Position in i-arm (0-based for structure)
                structure_pos = edit_pos - i_start
            # Check if edit is in j-arm
            elif j_start <= edit_pos <= j_end:
                # Position in j-arm, offset by i-arm length
                structure_pos = i_length + (edit_pos - j_start)
        
        else:  # strand == '-'
            # Reverse strand: coordinates go from high to low (3' to 5')
            # Need to reverse the mapping
            
            # For minus strand, genomic coords are reversed:
            # i_start is actually at the 3' end of i-arm
            # j_end is actually at the 5' end of j-arm
            
            # Check if edit is in i-arm (remember: i_end < i_start for minus)
            if i_end <= edit_pos <= i_start:
                # Position from 3' end of i-arm, but we need position from 5' end
                # So we reverse it
                pos_from_3prime = edit_pos - i_end
                structure_pos = i_length - 1 - pos_from_3prime
            
            # Check if edit is in j-arm (j_start < j_end for minus) 
            elif j_start <= edit_pos <= j_end:
                # Position from 3' end of j-arm
                pos_from_3prime = edit_pos - j_start
                # In structure, j-arm comes after i-arm
                # But for minus strand, we need to reverse within j-arm
                structure_pos = i_length + (j_length - 1 - pos_from_3prime)
        
        # Add to list if within structure bounds
        if structure_pos is not None and 0 <= structure_pos < (i_length + j_length):
            structure_edits.append([int(structure_pos), frequency])
    
    return structure_edits

def _safe_val(val, dtype=float, default=0):
    """Safely convert a pandas value, handling NaN"""
    if pd.isna(val):
        return default
    return dtype(val)


class DSRNARequestHandler(SimpleHTTPRequestHandler):
    """HTTP handler that serves the UI and API endpoints for browsing results"""

    # Class-level attributes set before server starts
    df = None
    editing_sites = None
    html_content = ""

    def end_headers(self):
        self.send_header('Access-Control-Allow-Origin', '*')
        super().end_headers()

    def log_message(self, format, *args):
        """Suppress default request logging"""
        pass

    def _send_json(self, data):
        body = json.dumps(data, ensure_ascii=False).encode('utf-8')
        self.send_response(200)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Content-Length', len(body))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        parsed = urlparse(self.path)
        path = parsed.path
        params = parse_qs(parsed.query)

        if path == '/' or path == '/index.html':
            body = self.html_content.encode('utf-8')
            self.send_response(200)
            self.send_header('Content-Type', 'text/html; charset=utf-8')
            self.send_header('Content-Length', len(body))
            self.end_headers()
            self.wfile.write(body)
        elif path == '/api/summary':
            self._handle_summary()
        elif path == '/api/results':
            self._handle_results(params)
        elif path.startswith('/api/result/'):
            try:
                result_id = int(path.split('/')[-1])
                self._handle_single_result(result_id)
            except (ValueError, IndexError):
                self.send_error(400, "Invalid result ID")
        elif path == '/api/download':
            self._handle_download(params)
        else:
            self.send_error(404)

    def _handle_summary(self):
        df = self.__class__.df
        chrom_counts = df.groupby('Chromosome').size().to_dict()
        # Natural sort chromosomes
        chroms = sorted(chrom_counts.keys(), key=_natural_sort_key)
        self._send_json({
            'total': len(df),
            'chromosomes': [{'name': c, 'count': chrom_counts[c]} for c in chroms],
            'has_editing': self.__class__.editing_sites is not None
        })

    def _get_filtered_df(self, params):
        """Apply all filters and return filtered DataFrame"""
        df = self.__class__.df

        # Search: supports "chr:start-end", chromosome name, or row index
        search = params.get('search', [None])[0]
        if search:
            import re
            # Try coordinate format: chr:start-end
            coord_match = re.match(r'^(.+?):(\d+)-(\d+)$', search.replace(',', ''))
            if coord_match:
                s_chr, s_start, s_end = coord_match.group(1), int(coord_match.group(2)), int(coord_match.group(3))
                df = df[(df['Chromosome'] == s_chr) &
                        (df['i_start'] <= s_end) & (df['j_end'] >= s_start)]
            elif search.isdigit():
                # Row index
                idx = int(search)
                if idx in df.index:
                    df = df.loc[[idx]]
                else:
                    df = df.iloc[0:0]  # empty
            else:
                # Chromosome name match (partial)
                df = df[df['Chromosome'].str.contains(search, case=False, na=False)]

        chrom = params.get('chromosome', [None])[0]
        strand = params.get('strand', [None])[0]
        if chrom:
            df = df[df['Chromosome'] == chrom]
        if strand:
            df = df[df['Strand'] == strand]
        if params.get('likely_edited', [None])[0] == '1' and 'likely_edited' in df.columns:
            df = df[df['likely_edited'].astype(str).str.strip().str.lower() == 'yes']
        if params.get('likely_forms', [None])[0] == '1' and 'likely_forms' in df.columns:
            df = df[df['likely_forms'].astype(str).str.strip().str.lower() == 'yes']

        # Column-level numeric filters (col_min, col_max)
        filter_map = {
            'score': 'Score', 'dG': 'dG(kcal/mol)', 'base_pairs': 'base_pairs',
            'percent_paired': 'percent_paired', 'longest_helix': 'longest_helix',
            'loop': None,
        }
        for key, col in filter_map.items():
            min_val = params.get(f'{key}_min', [None])[0]
            max_val = params.get(f'{key}_max', [None])[0]
            if key == 'loop':
                if min_val is not None or max_val is not None:
                    loop_col = df['j_start'] - df['i_end']
                    if min_val is not None:
                        df = df[loop_col >= float(min_val)]
                    if max_val is not None:
                        df = df[loop_col <= float(max_val)]
            elif col and col in df.columns:
                if min_val is not None:
                    df = df[df[col].astype(float) >= float(min_val)]
                if max_val is not None:
                    df = df[df[col].astype(float) <= float(max_val)]
        return df

    def _handle_results(self, params):
        df = self._get_filtered_df(params)
        total = len(df)

        # Sort
        sort_col = params.get('sort', ['Score'])[0]
        sort_dir = params.get('sort_dir', ['desc'])[0]
        col_map = {
            'score': 'Score', 'dG': 'dG(kcal/mol)', 'percent_paired': 'percent_paired',
            'longest_helix': 'longest_helix', 'i_start': 'i_start', 'base_pairs': 'base_pairs',
            'Score': 'Score', 'dG(kcal/mol)': 'dG(kcal/mol)',
            'stability_model_score': 'stability_model_score',
            'probing_model_score': 'probing_model_score'
        }
        actual_col = col_map.get(sort_col, 'Score')
        if actual_col in df.columns:
            df = df.sort_values(actual_col, ascending=(sort_dir == 'asc'), na_position='last')

        # Paginate
        page = int(params.get('page', [0])[0])
        page_size = int(params.get('page_size', [50])[0])
        page_size = min(page_size, 200)
        start = page * page_size
        page_df = df.iloc[start:start + page_size]

        # Return lightweight columns only (no sequences/structure)
        results = []
        for idx, row in page_df.iterrows():
            r = {
                'id': int(idx),
                'chromosome': str(row['Chromosome']),
                'strand': str(row['Strand']),
                'i_start': _safe_val(row['i_start'], int, 0),
                'i_end': _safe_val(row['i_end'], int, 0),
                'j_start': _safe_val(row['j_start'], int, 0),
                'j_end': _safe_val(row['j_end'], int, 0),
                'score': _safe_val(row['Score'], float, 0),
                'dG': _safe_val(row.get('dG(kcal/mol)', 0), float, 0),
                'percent_paired': _safe_val(row.get('percent_paired', 0), float, 0),
                'base_pairs': _safe_val(row.get('base_pairs', 0), int, 0),
                'longest_helix': _safe_val(row.get('longest_helix', 0), int, 0),
                'i_len': len(str(row['i_seq'])) if pd.notna(row.get('i_seq')) else 0,
                'j_len': len(str(row['j_seq'])) if pd.notna(row.get('j_seq')) else 0,
            }
            if 'stability_model_score' in row:
                r['stability_model_score'] = _safe_val(row['stability_model_score'], float, None)
            if 'probing_model_score' in row:
                r['probing_model_score'] = _safe_val(row['probing_model_score'], float, None)
            if 'likely_edited' in row:
                r['likely_edited'] = str(row['likely_edited']).strip().lower() == 'yes'
            if 'likely_forms' in row:
                r['likely_forms'] = str(row['likely_forms']).strip().lower() == 'yes'
            results.append(r)

        self._send_json({
            'total': total,
            'page': page,
            'page_size': page_size,
            'results': results
        })

    def _handle_download(self, params):
        """Download filtered results as TSV, applying the same filters as the table view"""
        # Reuse _handle_results filtering logic by extracting filtered df
        df = self._get_filtered_df(params)

        tsv = df.to_csv(sep='\t', index=False)
        body = tsv.encode('utf-8')
        self.send_response(200)
        self.send_header('Content-Type', 'text/tab-separated-values')
        self.send_header('Content-Disposition', 'attachment; filename="dsrnascan_filtered.tsv"')
        self.send_header('Content-Length', len(body))
        self.end_headers()
        self.wfile.write(body)

    def _handle_single_result(self, result_id):
        df = self.__class__.df
        if result_id not in df.index:
            self.send_error(404, "Result not found")
            return

        row = df.loc[result_id]
        result = {
            'id': result_id,
            'chromosome': str(row['Chromosome']),
            'strand': str(row['Strand']),
            'i_start': _safe_val(row['i_start'], int, 0),
            'i_end': _safe_val(row['i_end'], int, 0),
            'j_start': _safe_val(row['j_start'], int, 0),
            'j_end': _safe_val(row['j_end'], int, 0),
            'score': _safe_val(row['Score'], float, 0),
            'dG': _safe_val(row.get('dG(kcal/mol)', 0), float, 0),
            'percent_paired': _safe_val(row.get('percent_paired', 0), float, 0),
            'base_pairs': _safe_val(row.get('base_pairs', 0), int, 0),
            'longest_helix': _safe_val(row.get('longest_helix', 0), int, 0),
            'i_seq': str(row['i_seq']) if pd.notna(row['i_seq']) else '',
            'j_seq': str(row['j_seq']) if pd.notna(row['j_seq']) else '',
            'structure': str(row['structure']) if pd.notna(row['structure']) else '',
        }

        # Add effective coordinates
        for col in ['eff_i_start', 'eff_i_end', 'eff_j_start', 'eff_j_end']:
            if col in row and pd.notna(row[col]):
                result[col] = int(row[col])

        # Add ML scores and confidence labels
        for col in ['stability_model_score', 'probing_model_score']:
            if col in row and pd.notna(row[col]):
                result[col] = float(row[col])
        if 'likely_edited' in row:
            result['likely_edited'] = str(row['likely_edited']).strip().lower() == 'yes'
        if 'likely_forms' in row:
            result['likely_forms'] = str(row['likely_forms']).strip().lower() == 'yes'

        # Compute editing sites on demand
        editing_sites = self.__class__.editing_sites
        if editing_sites:
            result['editing_sites'] = map_editing_to_dsrna(row, editing_sites)

        self._send_json(result)


def _natural_sort_key(s):
    """Sort strings with embedded numbers naturally (chr1, chr2, ..., chr10)"""
    import re
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', str(s))]


def create_html_page(has_editing=False):
    """Generate the HTML page with API-driven UI"""
    return '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>dsRNA Structure Browser</title>
    <script src="https://d3js.org/d3.v3.min.js"></script>
    <script src="https://cdn.jsdelivr.net/gh/ViennaRNA/fornac@master/dist/fornac.js"></script>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/gh/ViennaRNA/fornac@master/dist/fornac.css">
    <style>
        * { box-sizing: border-box; margin: 0; padding: 0; }
        body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
               background: #f0f2f5; color: #1a1a2e; }
        .header { background: #1a1a2e; color: white; padding: 16px 24px; }
        .header h1 { font-size: 20px; font-weight: 600; }
        .header .subtitle { font-size: 13px; color: #8888aa; margin-top: 4px; }
        .controls { display: flex; gap: 12px; padding: 16px 24px; background: white;
                     border-bottom: 1px solid #e0e0e0; flex-wrap: wrap; align-items: center; }
        .controls select, .controls button { padding: 6px 12px; border: 1px solid #ccc;
                     border-radius: 4px; font-size: 13px; background: white; }
        .controls button { cursor: pointer; background: #f5f5f5; }
        .controls button:hover { background: #e8e8e8; }
        .controls .total { font-size: 13px; color: #666; margin-left: auto; }
        .main { display: flex; height: calc(100vh - 120px); }
        .table-panel { flex: 1; overflow-y: auto; border-right: 1px solid #e0e0e0; background: white; }
        .viewer-panel { flex: 1; display: flex; flex-direction: column; background: white; }
        table { width: 100%; border-collapse: collapse; font-size: 12px; }
        th { position: sticky; top: 0; background: #f8f9fa; padding: 8px 10px; text-align: left;
             border-bottom: 2px solid #dee2e6; cursor: pointer; white-space: nowrap; user-select: none; }
        th:hover { background: #e9ecef; }
        th.sorted-asc::after { content: " \\25B2"; font-size: 10px; }
        th.sorted-desc::after { content: " \\25BC"; font-size: 10px; }
        td { padding: 6px 10px; border-bottom: 1px solid #f0f0f0; white-space: nowrap; }
        tr { cursor: pointer; }
        tr:hover { background: #f0f4ff; }
        tr.selected { background: #d0e0ff; }
        .pagination { display: flex; align-items: center; justify-content: center; gap: 8px;
                       padding: 8px; border-top: 1px solid #e0e0e0; font-size: 13px; background: #fafafa; }
        .info-bar { padding: 12px 16px; background: #f8f9fa; border-bottom: 1px solid #e0e0e0;
                     font-size: 13px; display: none; }
        .info-bar .metrics { display: flex; gap: 24px; flex-wrap: wrap; }
        .info-bar .metric-label { color: #666; }
        .info-bar .metric-value { font-weight: 600; font-family: monospace; }
        #forna-viewer { flex: 1; min-height: 400px; }
        .loading { text-align: center; padding: 40px; color: #999; }
    </style>
</head>
<body>
    <div class="header">
        <h1>dsRNA Structure Browser</h1>
        <div class="subtitle" id="summary-text">Loading...</div>
    </div>
    <div class="controls">
        <input type="text" id="search-box" placeholder="Search: chr:start-end or ID" style="width:200px;padding:4px 8px;font-size:13px;">
        <label>Chr:</label>
        <select id="chrom-filter"><option value="">All</option></select>
        <label>Strand:</label>
        <select id="strand-filter">
            <option value="">Both</option>
            <option value="+">+</option>
            <option value="-">-</option>
        </select>
        <label>Sort:</label>
        <select id="sort-col">
            <option value="score">Score</option>
            <option value="dG">dG</option>
            <option value="percent_paired">% Paired</option>
            <option value="longest_helix">Helix</option>
            <option value="base_pairs">Base Pairs</option>
            <option value="i_start">Position</option>
            <option value="stability_model_score">Stability Model</option>
            <option value="probing_model_score">Probing Model</option>
        </select>
        <select id="sort-dir">
            <option value="desc">Desc</option>
            <option value="asc">Asc</option>
        </select>
        <label style="margin-left:12px;"><input type="checkbox" id="filter-edited"> Likely Edited</label>
        <label><input type="checkbox" id="filter-forms"> Likely Forms</label>
        <button id="download-btn" style="margin-left:12px;padding:4px 10px;cursor:pointer;">Download TSV</button>
        <div class="total" id="total-text"></div>
    </div>
    <div class="main">
        <div class="table-panel">
            <table>
                <thead><tr>
                    <th data-col="i_start">Location</th>
                    <th>Strand</th>
                    <th data-col="score">Score</th>
                    <th data-col="dG">dG</th>
                    <th data-col="base_pairs">BP</th>
                    <th data-col="percent_paired">%Paired</th>
                    <th data-col="longest_helix">Helix</th>
                    <th>Arms</th>
                    <th>Loop</th>
                    <th data-col="stability_model_score">Stability</th>
                    <th data-col="probing_model_score">Probing</th>
                </tr>
                <tr id="filter-row" style="background:#f0f0f0;">
                    <td></td><td></td>
                    <td><input type="number" class="col-filter" data-col="score" data-type="min" placeholder="min" style="width:40px;font-size:11px;"></td>
                    <td><input type="number" class="col-filter" data-col="dG" data-type="max" placeholder="max" style="width:45px;font-size:11px;" step="0.1"></td>
                    <td><input type="number" class="col-filter" data-col="base_pairs" data-type="min" placeholder="min" style="width:35px;font-size:11px;"></td>
                    <td><input type="number" class="col-filter" data-col="percent_paired" data-type="min" placeholder="min" style="width:35px;font-size:11px;"></td>
                    <td><input type="number" class="col-filter" data-col="longest_helix" data-type="min" placeholder="min" style="width:35px;font-size:11px;"></td>
                    <td></td>
                    <td><input type="number" class="col-filter" data-col="loop" data-type="max" placeholder="max" style="width:45px;font-size:11px;"></td>
                    <td></td><td></td>
                </tr>
                </thead>
                <tbody id="results-body"></tbody>
            </table>
            <div class="pagination">
                <button id="prev-btn" disabled>&lt; Prev</button>
                <span id="page-info">Page 1</span>
                <button id="next-btn">Next &gt;</button>
            </div>
        </div>
        <div class="viewer-panel">
            <div class="info-bar" id="info-bar">
                <div class="metrics" id="metrics"></div>
                <div id="seq-panel" style="display:none; margin-top:8px; padding-top:8px; border-top:1px solid #e0e0e0;"></div>
            </div>
            <div id="forna-viewer"><div class="loading">Select a structure from the table</div></div>
        </div>
    </div>
<script>
let state = { page: 0, pageSize: 50, chromosome: '', strand: '', sort: 'score', sortDir: 'desc', total: 0, likelyEdited: false, likelyForms: false, filters: {}, search: '' };

async function fetchJSON(url) { return (await fetch(url)).json(); }

async function init() {
    const summary = await fetchJSON('/api/summary');
    document.getElementById('summary-text').textContent =
        `${summary.total.toLocaleString()} dsRNA structures across ${summary.chromosomes.length} chromosomes`;
    const sel = document.getElementById('chrom-filter');
    summary.chromosomes.forEach(c => {
        const opt = document.createElement('option');
        opt.value = c.name;
        opt.textContent = `${c.name} (${c.count.toLocaleString()})`;
        sel.appendChild(opt);
    });
    loadResults();
}

async function loadResults() {
    const params = new URLSearchParams({
        page: state.page, page_size: state.pageSize,
        sort: state.sort, sort_dir: state.sortDir
    });
    if (state.chromosome) params.set('chromosome', state.chromosome);
    if (state.strand) params.set('strand', state.strand);
    if (state.likelyEdited) params.set('likely_edited', '1');
    if (state.likelyForms) params.set('likely_forms', '1');
    if (state.search) params.set('search', state.search);
    for (const [key, val] of Object.entries(state.filters)) {
        params.set(key, val);
    }

    const data = await fetchJSON('/api/results?' + params);
    state.total = data.total;

    const totalPages = Math.ceil(data.total / state.pageSize);
    document.getElementById('total-text').textContent = `${data.total.toLocaleString()} results`;
    document.getElementById('page-info').textContent = `Page ${state.page + 1} of ${totalPages || 1}`;
    document.getElementById('prev-btn').disabled = state.page === 0;
    document.getElementById('next-btn').disabled = state.page >= totalPages - 1;

    const tbody = document.getElementById('results-body');
    tbody.innerHTML = '';
    data.results.forEach(r => {
        const tr = document.createElement('tr');
        tr.dataset.id = r.id;
        const loop = r.j_start - r.i_end;
        function scoreColor(val, threshold) {
            if (val == null) return 'color:#ccc';
            if (val >= threshold) {
                const t = Math.min((val - threshold) / (1 - threshold), 1);
                const g = Math.round(100 + t * 80);
                return 'color:rgb(0,' + g + ',0);font-weight:bold';
            } else {
                return 'color:#ccc';
            }
        }
        const stabVal = r.stability_model_score != null ? r.stability_model_score.toFixed(3) : '-';
        const probVal = r.probing_model_score != null ? r.probing_model_score.toFixed(3) : '-';
        const stabStyle = scoreColor(r.stability_model_score, 0.247);
        const probStyle = scoreColor(r.probing_model_score, 0.032);
        tr.innerHTML = `<td>${r.chromosome}:${r.i_start.toLocaleString()}-${r.j_end.toLocaleString()}</td>
            <td>${r.strand}</td>
            <td>${r.score}</td><td>${r.dG.toFixed(1)}</td><td>${r.base_pairs}</td>
            <td>${r.percent_paired.toFixed(1)}%</td><td>${r.longest_helix}</td>
            <td>${r.i_len}/${r.j_len}</td><td>${loop.toLocaleString()}</td>
            <td style="${stabStyle}">${stabVal}</td><td style="${probStyle}">${probVal}</td>`;
        tr.addEventListener('click', () => selectResult(r.id, tr));
        tbody.appendChild(tr);
    });
}

async function selectResult(id, tr) {
    document.querySelectorAll('tr.selected').forEach(el => el.classList.remove('selected'));
    if (tr) tr.classList.add('selected');

    const data = await fetchJSON('/api/result/' + id);
    showInfo(data);
    showStructure(data);
}

function showInfo(d) {
    const bar = document.getElementById('info-bar');
    bar.style.display = 'block';
    let html = `<div><span class="metric-label">Location:</span> <span class="metric-value">${d.chromosome}:${d.i_start.toLocaleString()}-${d.j_end.toLocaleString()} (${d.strand})</span></div>
        <div><span class="metric-label">Energy:</span> <span class="metric-value">${d.dG.toFixed(1)} kcal/mol</span></div>
        <div><span class="metric-label">Score:</span> <span class="metric-value">${d.score}</span></div>
        <div><span class="metric-label">Base Pairs:</span> <span class="metric-value">${d.base_pairs} (${d.percent_paired.toFixed(1)}%)</span></div>
        <div><span class="metric-label">Helix:</span> <span class="metric-value">${d.longest_helix} bp</span></div>`;
    if (d.stability_model_score != null) {
        if (d.likely_edited) {
            html += '<div><span class="metric-label">Stability Model:</span> <span class="metric-value">' + d.stability_model_score.toFixed(3) + '</span> <span style="color:#27ae60;font-size:12px;font-weight:bold;">high confidence editing</span></div>';
        } else {
            html += '<div><span class="metric-label">Stability Model:</span> <span class="metric-value">' + d.stability_model_score.toFixed(3) + '</span> <span style="color:#e74c3c;font-size:12px;">low confidence editing</span></div>';
        }
    }
    if (d.probing_model_score != null) {
        if (d.likely_forms) {
            html += '<div><span class="metric-label">Probing Model:</span> <span class="metric-value">' + d.probing_model_score.toFixed(3) + '</span> <span style="color:#2980b9;font-size:12px;font-weight:bold;">high confidence formation</span></div>';
        } else {
            html += '<div><span class="metric-label">Probing Model:</span> <span class="metric-value">' + d.probing_model_score.toFixed(3) + '</span> <span style="color:#e74c3c;font-size:12px;">low confidence formation</span></div>';
        }
    }
    html += '<div style="margin-top:6px;font-size:10px;color:#888;line-height:1.4;">'
        + 'Stability model predicts ADAR editing substrate likelihood from structural features (threshold: 0.247). '
        + 'Probing model predicts in vivo structure formation from RNA probing data (threshold: 0.032). '
        + 'Thresholds derived from human editing sites and may not apply to other organisms.</div>';
    if (d.editing_sites && d.editing_sites.length > 0) html += `<div><span class="metric-label" style="color:#27AE60">Editing Sites:</span> <span class="metric-value">${d.editing_sites.length}</span></div>`;
    document.getElementById('metrics').innerHTML = html;

    // Populate sequence panel - hidden data with copy buttons only
    const seqPanel = document.getElementById('seq-panel');
    if (d.i_seq && d.j_seq) {
        const parts = d.structure ? d.structure.split('&') : ['',''];
        const btnStyle = 'font-size:11px; padding:2px 8px; cursor:pointer; margin-left:6px; border:1px solid #ccc; border-radius:3px; background:#fff;';
        const dbnHeader = `>${d.id} ${d.location} dG=${d.energy}`;
        const dbnSeq = d.i_seq + 'N' + d.j_seq;
        const dbnStr = parts[0] + '.' + parts[1];
        const dbnFull = dbnHeader + '\\n' + dbnSeq + '\\n' + dbnStr;
        seqPanel.style.display = 'block';
        seqPanel.innerHTML = `
            <span style="font-family:sans-serif; font-size:12px; color:#666;"><strong>Copy:</strong></span>
            <button onclick="copyText('dbn-full',this)" style="${btnStyle}">DBN</button>
            <button onclick="copyText('seq-i',this)" style="${btnStyle}">i-seq (${d.i_seq.length} nt)</button>
            <button onclick="copyText('str-i',this)" style="${btnStyle}">i-struct</button>
            <button onclick="copyText('seq-j',this)" style="${btnStyle}">j-seq (${d.j_seq.length} nt)</button>
            <button onclick="copyText('str-j',this)" style="${btnStyle}">j-struct</button>
            <span id="dbn-full" style="display:none;">${dbnFull}</span>
            <span id="seq-i" style="display:none;">${d.i_seq}</span>
            <span id="str-i" style="display:none;">${parts[0]}</span>
            <span id="seq-j" style="display:none;">${d.j_seq}</span>
            <span id="str-j" style="display:none;">${parts[1]}</span>`;
    } else {
        seqPanel.style.display = 'none';
    }
}

function showStructure(data) {
    const viewer = document.getElementById('forna-viewer');
    viewer.innerHTML = '';

    const parts = data.structure.split('&');
    const sequence = data.i_seq + data.j_seq;
    const structure = parts[0] + parts[1];
    const maxArm = Math.max(data.i_seq.length, data.j_seq.length);

    // For large structures, skip Forna rendering (sequences available via copy buttons)
    if (maxArm > 500) {
        viewer.innerHTML = `<div class="loading">Arms too large for interactive view (${maxArm} nt). Use copy buttons above.</div>`;
        return;
    }

    const container = document.createElement('div');
    container.id = 'forna-container';
    container.style.width = '100%';
    container.style.height = '100%';
    viewer.appendChild(container);

    const rnaViz = new fornac.FornaContainer("#forna-container", {
        applyForce: true, allowPanningAndZooming: true,
        initialSize: [viewer.offsetWidth || 600, viewer.offsetHeight || 500],
        friction: 0.35, middleCharge: -30, otherCharge: -30
    });
    rnaViz.addRNA(structure, { structure: structure, sequence: sequence, labelInterval: 10 });

    if (data.editing_sites && data.editing_sites.length > 0) {
        setTimeout(() => annotateEditingSites('forna-container', data.editing_sites), 1500);
    }
    setTimeout(() => { rnaViz.centerView(); rnaViz.zoomToFit(); }, 1000);
}

function copyText(elemId, btn) {
    const text = document.getElementById(elemId).textContent;
    navigator.clipboard.writeText(text).then(() => {
        const orig = btn.textContent;
        btn.textContent = 'Copied!';
        setTimeout(() => { btn.textContent = orig; }, 1500);
    });
}

function annotateEditingSites(containerId, editingSites) {
    const svg = d3.select('#' + containerId + ' svg');
    const allCircles = svg.selectAll('circle');
    const radiusCounts = {};
    allCircles.each(function() { const r = parseFloat(d3.select(this).attr('r')); radiusCounts[r] = (radiusCounts[r]||0)+1; });
    let mostCommonRadius = 0, maxCount = 0;
    for (const [r, c] of Object.entries(radiusCounts)) { if (c > maxCount) { maxCount = c; mostCommonRadius = parseFloat(r); } }
    const nucs = [];
    allCircles.each(function() { if (Math.abs(parseFloat(d3.select(this).attr('r')) - mostCommonRadius) < 0.1) nucs.push(this); });
    editingSites.forEach(site => {
        const [pos, freq] = Array.isArray(site) ? site : [site, 1.0];
        if (pos < nucs.length) {
            const circle = d3.select(nucs[pos]);
            const color = freq >= 0.8 ? '#0B5345' : freq >= 0.5 ? '#148F77' : freq >= 0.3 ? '#27AE60' : '#52BE80';
            circle.style('stroke', color).style('stroke-width', '3px').style('stroke-opacity', 1);
            let title = circle.select('title');
            if (title.empty()) title = circle.append('title');
            title.text('Editing site pos ' + (pos+1) + ' (' + (freq*100).toFixed(1) + '%)');
        }
    });
}

// Column filter event listeners
document.querySelectorAll('.col-filter').forEach(input => {
    input.addEventListener('change', () => {
        state.filters = {};
        document.querySelectorAll('.col-filter').forEach(el => {
            if (el.value !== '') {
                const key = el.dataset.col + '_' + el.dataset.type;
                state.filters[key] = parseFloat(el.value);
            }
        });
        state.page = 0;
        loadResults();
    });
});

// Search box
document.getElementById('search-box').addEventListener('keydown', e => {
    if (e.key === 'Enter') {
        state.search = e.target.value.trim();
        state.page = 0;
        loadResults();
    }
});

// Event listeners
document.getElementById('chrom-filter').addEventListener('change', e => { state.chromosome = e.target.value; state.page = 0; loadResults(); });
document.getElementById('strand-filter').addEventListener('change', e => { state.strand = e.target.value; state.page = 0; loadResults(); });
document.getElementById('sort-col').addEventListener('change', e => { state.sort = e.target.value; state.page = 0; loadResults(); });
document.getElementById('sort-dir').addEventListener('change', e => { state.sortDir = e.target.value; state.page = 0; loadResults(); });
document.getElementById('filter-edited').addEventListener('change', e => { state.likelyEdited = e.target.checked; state.page = 0; loadResults(); });
document.getElementById('filter-forms').addEventListener('change', e => { state.likelyForms = e.target.checked; state.page = 0; loadResults(); });
document.getElementById('prev-btn').addEventListener('click', () => { if (state.page > 0) { state.page--; loadResults(); } });
document.getElementById('next-btn').addEventListener('click', () => { state.page++; loadResults(); });
document.getElementById('download-btn').addEventListener('click', () => {
    const params = new URLSearchParams();
    if (state.chromosome) params.set('chromosome', state.chromosome);
    if (state.strand) params.set('strand', state.strand);
    if (state.likelyEdited) params.set('likely_edited', '1');
    if (state.likelyForms) params.set('likely_forms', '1');
    for (const [key, val] of Object.entries(state.filters)) {
        params.set(key, val);
    }
    window.location.href = '/api/download?' + params;
});

// Column header sorting
document.querySelectorAll('th[data-col]').forEach(th => {
    th.addEventListener('click', () => {
        const col = th.dataset.col;
        if (state.sort === col) { state.sortDir = state.sortDir === 'desc' ? 'asc' : 'desc'; }
        else { state.sort = col; state.sortDir = 'desc'; }
        document.getElementById('sort-col').value = state.sort;
        document.getElementById('sort-dir').value = state.sortDir;
        state.page = 0;
        loadResults();
    });
});

init();
</script>
</body>
</html>'''


def main():
    parser = argparse.ArgumentParser(description='Browse dsRNAscan results with Forna visualization')
    parser.add_argument('output_dir', nargs='?', default='.',
                       help='dsRNAscan output directory (default: current directory)')
    parser.add_argument('--editing-file', type=str,
                       help='BED or GFF3 file with RNA editing sites')
    parser.add_argument('--large-editing-file', action='store_true',
                       help='Force optimized parsing for large editing files')
    parser.add_argument('--port', type=int, default=8080,
                       help='Port for the web server (default: 8080)')
    parser.add_argument('--no-browser', action='store_true',
                       help='Do not automatically open the web browser')

    args = parser.parse_args()

    # Find results files
    output_path = Path(args.output_dir)
    if not output_path.exists():
        print(f"Error: Directory {args.output_dir} does not exist")
        sys.exit(1)

    results_files = list(output_path.glob('*_merged_results.txt'))
    if not results_files:
        print(f"Error: No dsRNAscan results files (*_merged_results.txt) found in {args.output_dir}")
        sys.exit(1)

    # Load all results into a single DataFrame
    dfs = []
    for rf in results_files:
        print(f"Loading {rf}...")
        try:
            df = pd.read_csv(rf, sep='\t')
            df = df.dropna(subset=['Chromosome'])
            dfs.append(df)
        except Exception as e:
            print(f"  Warning: {e}")

    if not dfs:
        print("Error: No valid results found")
        sys.exit(1)

    all_df = pd.concat(dfs, ignore_index=True)
    print(f"Loaded {len(all_df):,} dsRNA structures")

    # Parse editing sites if provided
    editing_sites = None
    if args.editing_file:
        needed_chroms = set(all_df['Chromosome'].unique())
        print(f"Loading editing sites from {args.editing_file}...")
        file_size = os.path.getsize(args.editing_file)
        use_optimized = args.large_editing_file or file_size > 10_000_000
        if use_optimized:
            editing_sites = parse_editing_file_optimized(args.editing_file, needed_chroms)
        else:
            editing_sites = parse_editing_file_optimized(args.editing_file, None)

    # Set up handler with data
    DSRNARequestHandler.df = all_df
    DSRNARequestHandler.editing_sites = editing_sites
    DSRNARequestHandler.html_content = create_html_page(has_editing=bool(editing_sites))

    # Start server
    with socketserver.TCPServer(("", args.port), DSRNARequestHandler) as httpd:
        print(f"\nBrowser: http://localhost:{args.port}/")
        print("Press Ctrl+C to stop")

        if not args.no_browser:
            webbrowser.open(f'http://localhost:{args.port}/')

        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nStopped.")


if __name__ == '__main__':
    main()