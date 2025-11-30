#!/usr/bin/env python3
"""
Consolidate protein design metrics from ipSAE, Prodigy, Foldseek, and sequence outputs.

Generates an interactive HTML report with AG Grid for searching/filtering.

This script reads from staged files in a flat directory structure.
"""

import argparse
import csv
import os
import sys
import re
import json
from pathlib import Path
from collections import defaultdict


def parse_ipsae_file(ipsae_file):
    """
    Parse ipSAE output file and extract the max ipSAE score.

    Returns the ipSAE value from the row where Type='max'.
    """
    if not os.path.exists(ipsae_file):
        return None

    try:
        with open(ipsae_file, 'r') as f:
            lines = f.readlines()

        # Find header line to get column indices
        header_idx = None
        for i, line in enumerate(lines):
            if 'ipSAE' in line and 'Type' in line:
                header_idx = i
                break

        if header_idx is None:
            return None

        # Parse header
        header = lines[header_idx].split()
        try:
            type_col = header.index('Type')
            ipsae_col = header.index('ipSAE')
        except ValueError:
            return None

        # Find max row
        for line in lines[header_idx + 1:]:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) > max(type_col, ipsae_col):
                if parts[type_col] == 'max':
                    return float(parts[ipsae_col])

    except Exception as e:
        print(f"Warning: Could not parse ipSAE file {ipsae_file}: {e}", file=sys.stderr)

    return None


def parse_prodigy_results(prodigy_file):
    """
    Parse Prodigy results.txt and extract the predicted binding affinity.

    Returns the value from "Predicted binding affinity (kcal.mol-1):" line.
    """
    if not os.path.exists(prodigy_file):
        return None

    try:
        with open(prodigy_file, 'r') as f:
            for line in f:
                if 'Predicted binding affinity (kcal.mol-1):' in line:
                    # Extract the value after the colon
                    match = re.search(r'Predicted binding affinity \(kcal\.mol-1\):\s*([-\d.]+)', line)
                    if match:
                        return float(match.group(1))
    except Exception as e:
        print(f"Warning: Could not parse Prodigy file {prodigy_file}: {e}", file=sys.stderr)

    return None


def parse_foldseek_summary(foldseek_file):
    """
    Parse Foldseek summary TSV and find the hit with shortest AA distance.

    Distance formula: Differences = Alignment Length × (1.0 - fident)

    Returns tuple: (target_name, distance) or (None, None)
    """
    if not os.path.exists(foldseek_file):
        return None, None

    try:
        best_target = None
        best_distance = float('inf')

        with open(foldseek_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) < 4:
                    continue

                # Column indices (0-based):
                # 0: query, 1: target, 2: fident, 3: alnlen
                try:
                    target = parts[1]
                    fident = float(parts[2])
                    alnlen = int(parts[3])

                    # Calculate distance (AA differences)
                    distance = round(alnlen * (1.0 - fident))

                    if distance < best_distance:
                        best_distance = distance
                        best_target = target
                except (ValueError, IndexError):
                    continue

        if best_target:
            return best_target, int(best_distance)

    except Exception as e:
        print(f"Warning: Could not parse Foldseek file {foldseek_file}: {e}", file=sys.stderr)

    return None, None


def parse_afdb_id(foldseek_target):
    """
    Parse AlphaFold DB ID from Foldseek target name.

    Input formats:
        - AF-P12345-F1-model_v4
        - AF-Q9Y6K9-F1-model_v4

    Returns tuple: (uniprot_id, afdb_url) or (None, None)
    """
    if not foldseek_target:
        return None, None

    # Pattern: AF-{UniProt}-F{fragment}-model_v{version}
    match = re.match(r'^AF-([A-Z0-9]+)-F\d+-model_v\d+$', foldseek_target)
    if match:
        uniprot_id = match.group(1)
        afdb_url = f"https://alphafold.ebi.ac.uk/entry/{uniprot_id}"
        return uniprot_id, afdb_url

    return None, None


def parse_sequence_file(seq_file):
    """
    Parse FASTA sequence file and extract the binder sequence.

    ProteinMPNN sequences may contain both binder and target separated by '/'.
    This function extracts only the binder sequence (before the '/').

    Returns the amino acid sequence string or None.
    """
    if not os.path.exists(seq_file):
        return None

    try:
        sequences = []
        current_header = None
        current_seq = []

        with open(seq_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_header and current_seq:
                        sequences.append((current_header, ''.join(current_seq)))
                    current_header = line
                    current_seq = []
                else:
                    current_seq.append(line)

            if current_header and current_seq:
                sequences.append((current_header, ''.join(current_seq)))

        if not sequences:
            return None

        # Get the shortest sequence (the binder is the designed smaller protein)
        full_sequence = min(sequences, key=lambda x: len(x[1]))[1]
        
        # Extract only the binder sequence (before the '/' if present)
        # ProteinMPNN format: BINDER_SEQUENCE/TARGET_SEQUENCE
        binder_sequence = full_sequence.split('/')[0]
        
        return binder_sequence

    except Exception as e:
        print(f"Warning: Could not parse sequence file {seq_file}: {e}", file=sys.stderr)

    return None


def extract_design_id_from_ipsae(filename):
    """
    Extract design ID from ipSAE filename.

    Boltz2 CIF files are named: {design_id}_model_0.cif
    ipSAE outputs: {design_id}_model_0_{pae}_{dist}.txt

    Example: 2vsm_r1_s0_model_0_10_10.txt -> 2vsm_r1_s0
    """
    base = Path(filename).stem
    # Remove _XX_XX cutoff pattern at the end (e.g., _10_10)
    base = re.sub(r'_\d+_\d+$', '', base)
    # Remove _model_X pattern (from Boltz2 output naming)
    base = re.sub(r'_model_\d+$', '', base)
    return base


def extract_design_id_from_prodigy(filename):
    """
    Extract design ID from Prodigy filename.

    Example: 2vsm_r1_s0_prodigy_results.txt -> 2vsm_r1_s0
    """
    base = Path(filename).stem
    return base.replace('_prodigy_results', '')


def extract_design_id_from_foldseek(filename):
    """
    Extract design ID from Foldseek filename.

    Example: 2vsm_r1_s0_foldseek_summary.tsv -> 2vsm_r1_s0
    """
    base = Path(filename).stem
    return base.replace('_foldseek_summary', '').replace('_foldseek_results', '')


def extract_design_id_from_sequence(filename):
    """
    Extract design ID from sequence filename.

    PREPARE_BOLTZ2_SEQUENCES outputs: {sample}_r{rank}_s{seq_num}.fa
    The filename stem IS the design ID (no conversion needed).

    Example: 2vsm_r1_s0.fa -> 2vsm_r1_s0
    """
    return Path(filename).stem


def collect_metrics_from_dirs(ipsae_dir, prodigy_dir, foldseek_dir, sequence_dir=None, ipsae_cutoffs="10_10"):
    """
    Collect metrics from staged files in separate directories.

    Args:
        ipsae_dir: Directory containing ipSAE score files
        prodigy_dir: Directory containing Prodigy result files
        foldseek_dir: Directory containing Foldseek summary files
        sequence_dir: Directory containing sequence FASTA files
        ipsae_cutoffs: Suffix for ipSAE files (e.g., "10_10")

    Returns a dict mapping design_id -> metrics dict
    """
    metrics = defaultdict(dict)

    print(f"Using ipSAE cutoffs: {ipsae_cutoffs}")

    # 1. Parse ipSAE files (pattern: *_XX_XX.txt, excluding *_byres.txt)
    ipsae_path = Path(ipsae_dir)
    if ipsae_path.exists():
        ipsae_pattern = f'*_{ipsae_cutoffs}.txt'
        ipsae_files = list(ipsae_path.glob(ipsae_pattern))
        print(f"Found {len(ipsae_files)} ipSAE files in {ipsae_dir}")

        for ipsae_file in ipsae_files:
            if 'byres' in ipsae_file.name:
                continue

            design_id = extract_design_id_from_ipsae(ipsae_file.name)
            ipsae_score = parse_ipsae_file(str(ipsae_file))

            if ipsae_score is not None:
                metrics[design_id]['ipsae'] = ipsae_score
                print(f"  ipSAE [{design_id}]: {ipsae_score}")
    else:
        print(f"ipSAE directory not found: {ipsae_dir}")

    # 2. Parse Prodigy files (pattern: *_prodigy_results.txt)
    prodigy_path = Path(prodigy_dir)
    if prodigy_path.exists():
        prodigy_files = list(prodigy_path.glob('*_prodigy_results.txt'))
        print(f"Found {len(prodigy_files)} Prodigy files in {prodigy_dir}")

        for prodigy_file in prodigy_files:
            design_id = extract_design_id_from_prodigy(prodigy_file.name)
            binding_affinity = parse_prodigy_results(str(prodigy_file))

            if binding_affinity is not None:
                metrics[design_id]['binding_affinity'] = binding_affinity
                print(f"  Prodigy [{design_id}]: {binding_affinity}")
    else:
        print(f"Prodigy directory not found: {prodigy_dir}")

    # 3. Parse Foldseek files (pattern: *_foldseek_summary.tsv)
    foldseek_path = Path(foldseek_dir)
    if foldseek_path.exists():
        foldseek_files = list(foldseek_path.glob('*_foldseek_summary.tsv'))
        print(f"Found {len(foldseek_files)} Foldseek files in {foldseek_dir}")

        for foldseek_file in foldseek_files:
            design_id = extract_design_id_from_foldseek(foldseek_file.name)
            target, distance = parse_foldseek_summary(str(foldseek_file))

            if target is not None:
                metrics[design_id]['foldseek_target'] = target
                metrics[design_id]['foldseek_distance'] = distance
                # Parse AFDB ID
                uniprot_id, afdb_url = parse_afdb_id(target)
                if uniprot_id:
                    metrics[design_id]['afdb_id'] = uniprot_id
                    metrics[design_id]['afdb_url'] = afdb_url
                print(f"  Foldseek [{design_id}]: {target} (dist={distance}, uniprot={uniprot_id})")
    else:
        print(f"Foldseek directory not found: {foldseek_dir}")

    # 4. Parse sequence files (pattern: *.fa or *.fasta)
    if sequence_dir:
        sequence_path = Path(sequence_dir)
        if sequence_path.exists():
            seq_files = list(sequence_path.glob('*.fa')) + list(sequence_path.glob('*.fasta'))
            print(f"Found {len(seq_files)} sequence files in {sequence_dir}")

            for seq_file in seq_files:
                design_id = extract_design_id_from_sequence(seq_file.name)
                sequence = parse_sequence_file(str(seq_file))

                if sequence:
                    metrics[design_id]['sequence'] = sequence
                    metrics[design_id]['aa_length'] = len(sequence)
                    print(f"  Sequence [{design_id}]: {len(sequence)} AA")
        else:
            print(f"Sequence directory not found: {sequence_dir}")

    # ========================================
    # Summary: Show ID matching status
    # ========================================
    print("\n" + "=" * 60)
    print("DESIGN ID MATCHING SUMMARY")
    print("=" * 60)
    for design_id in sorted(metrics.keys()):
        data = metrics[design_id]
        has_ipsae = 'ipsae' in data
        has_prodigy = 'binding_affinity' in data
        has_foldseek = 'foldseek_target' in data
        has_sequence = 'sequence' in data

        status = []
        if has_ipsae: status.append("ipSAE")
        if has_prodigy: status.append("Prodigy")
        if has_foldseek: status.append("Foldseek")
        if has_sequence: status.append("Sequence")

        missing = []
        if not has_ipsae: missing.append("ipSAE")
        if not has_prodigy: missing.append("Prodigy")
        if not has_foldseek: missing.append("Foldseek")
        if not has_sequence: missing.append("Sequence")

        status_str = ", ".join(status) if status else "NONE"
        missing_str = f" [MISSING: {', '.join(missing)}]" if missing else " [COMPLETE]"

        print(f"  {design_id}: {status_str}{missing_str}")

    print("=" * 60 + "\n")

    return metrics


def generate_html_report(metrics, output_file, title="Protein Design Metrics Report"):
    """
    Generate an interactive HTML report with AG Grid for searching/filtering.
    """
    # Convert metrics dict to sorted list
    rows = []
    for design_id, data in sorted(metrics.items()):
        row = {
            'design_id': design_id,
            'ipsae': data.get('ipsae'),
            'binding_affinity': data.get('binding_affinity'),
            'foldseek_target': data.get('foldseek_target'),
            'foldseek_distance': data.get('foldseek_distance'),
            'afdb_id': data.get('afdb_id'),
            'afdb_url': data.get('afdb_url'),
            'sequence': data.get('sequence'),
            'aa_length': data.get('aa_length')
        }
        rows.append(row)

    # Sort by ipSAE ascending (lower = better interface prediction)
    rows.sort(key=lambda x: x.get('ipsae') or float('inf'))

    # Convert rows to JSON for AG Grid
    rows_json = json.dumps(rows)

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>

    <!-- AG Grid CSS -->
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/ag-grid-community@32.3.3/styles/ag-grid.css">
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/ag-grid-community@32.3.3/styles/ag-theme-alpine.css">

    <style>
        * {{
            box-sizing: border-box;
        }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f0f2f5;
        }}

        .container {{
            max-width: 1600px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.08);
        }}

        h1 {{
            color: #1a1a2e;
            margin-bottom: 8px;
            font-size: 1.8rem;
            font-weight: 600;
        }}

        .summary {{
            color: #666;
            margin-bottom: 24px;
            padding: 16px 20px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            border-radius: 8px;
            color: white;
            display: flex;
            gap: 30px;
            flex-wrap: wrap;
        }}

        .summary-item {{
            display: flex;
            flex-direction: column;
        }}

        .summary-label {{
            font-size: 0.85rem;
            opacity: 0.9;
        }}

        .summary-value {{
            font-size: 1.5rem;
            font-weight: 600;
        }}

        .metric-info {{
            font-size: 0.9em;
            color: #666;
            margin-bottom: 20px;
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
            gap: 12px;
        }}

        .metric-card {{
            background: #f8f9fa;
            padding: 12px 16px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}

        .metric-card dt {{
            font-weight: 600;
            color: #1a1a2e;
            margin-bottom: 4px;
        }}

        .metric-card dd {{
            margin: 0;
            color: #666;
            font-size: 0.85rem;
        }}

        .toolbar {{
            display: flex;
            gap: 12px;
            margin-bottom: 16px;
            flex-wrap: wrap;
            align-items: center;
        }}

        .search-box {{
            flex: 1;
            min-width: 200px;
            max-width: 400px;
        }}

        .search-box input {{
            width: 100%;
            padding: 10px 16px;
            border: 2px solid #e0e0e0;
            border-radius: 8px;
            font-size: 0.95rem;
            transition: border-color 0.2s;
        }}

        .search-box input:focus {{
            outline: none;
            border-color: #667eea;
        }}

        .btn {{
            padding: 10px 20px;
            border: none;
            border-radius: 8px;
            font-size: 0.9rem;
            font-weight: 500;
            cursor: pointer;
            transition: all 0.2s;
        }}

        .btn-primary {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
        }}

        .btn-primary:hover {{
            transform: translateY(-1px);
            box-shadow: 0 4px 12px rgba(102, 126, 234, 0.4);
        }}

        .btn-secondary {{
            background: #f0f0f0;
            color: #333;
        }}

        .btn-secondary:hover {{
            background: #e0e0e0;
        }}

        #metricsGrid {{
            height: 600px;
            width: 100%;
        }}

        .ag-theme-alpine {{
            --ag-header-background-color: #667eea;
            --ag-header-foreground-color: white;
            --ag-row-hover-color: #f0f4ff;
            --ag-selected-row-background-color: #e8f0fe;
            --ag-font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
        }}

        .ag-theme-alpine .ag-header-cell {{
            font-weight: 600;
        }}

        .good-value {{
            color: #28a745;
            font-weight: 600;
        }}

        .warning-value {{
            color: #ffc107;
        }}

        .afdb-link {{
            color: #667eea;
            text-decoration: none;
            font-weight: 500;
        }}

        .afdb-link:hover {{
            text-decoration: underline;
        }}

        .sequence-cell {{
            font-family: 'Monaco', 'Menlo', 'Consolas', monospace;
            font-size: 0.8rem;
            cursor: pointer;
        }}

        .sequence-collapsed {{
            max-width: 120px;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
            background: #f0f0f0;
            padding: 4px 8px;
            border-radius: 4px;
            display: inline-block;
        }}

        .sequence-collapsed:hover {{
            background: #e0e0e0;
        }}

        .sequence-expanded {{
            word-break: break-all;
            white-space: pre-wrap;
            background: #f8f9fa;
            padding: 8px;
            border-radius: 4px;
            max-height: 200px;
            overflow-y: auto;
            border: 1px solid #e0e0e0;
        }}

        .expand-icon {{
            margin-left: 4px;
            font-size: 0.7rem;
            opacity: 0.6;
        }}

        /* Modal for sequence */
        .modal {{
            display: none;
            position: fixed;
            z-index: 1000;
            left: 0;
            top: 0;
            width: 100%;
            height: 100%;
            background-color: rgba(0,0,0,0.5);
        }}

        .modal-content {{
            background-color: white;
            margin: 10% auto;
            padding: 24px;
            border-radius: 12px;
            width: 80%;
            max-width: 800px;
            max-height: 70vh;
            overflow-y: auto;
        }}

        .modal-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 16px;
        }}

        .modal-header h3 {{
            margin: 0;
            color: #1a1a2e;
        }}

        .close-btn {{
            background: none;
            border: none;
            font-size: 1.5rem;
            cursor: pointer;
            color: #666;
        }}

        .close-btn:hover {{
            color: #333;
        }}

        .sequence-display {{
            font-family: 'Monaco', 'Menlo', 'Consolas', monospace;
            font-size: 0.9rem;
            word-break: break-all;
            white-space: pre-wrap;
            background: #f8f9fa;
            padding: 16px;
            border-radius: 8px;
            border: 1px solid #e0e0e0;
            line-height: 1.6;
        }}

        .copy-btn {{
            margin-top: 12px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>{title}</h1>

        <div class="summary">
            <div class="summary-item">
                <span class="summary-label">Total Designs</span>
                <span class="summary-value">{len(rows)}</span>
            </div>
            <div class="summary-item">
                <span class="summary-label">With Sequences</span>
                <span class="summary-value">{sum(1 for r in rows if r.get('sequence'))}</span>
            </div>
            <div class="summary-item">
                <span class="summary-label">AFDB Matches</span>
                <span class="summary-value">{sum(1 for r in rows if r.get('afdb_id'))}</span>
            </div>
        </div>

        <div class="metric-info">
            <div class="metric-card">
                <dl>
                    <dt>ipSAE</dt>
                    <dd>Interface interaction prediction Score from Aligned Errors. Higher is better.</dd>
                </dl>
            </div>
            <div class="metric-card">
                <dl>
                    <dt>Binding Affinity</dt>
                    <dd>Predicted from PRODIGY (kcal/mol). More negative is better; &lt;-10 indicates strong binding.</dd>
                </dl>
            </div>
            <div class="metric-card">
                <dl>
                    <dt>Foldseek Distance</dt>
                    <dd>AA differences to closest structural match. Lower = more similar to known structures.</dd>
                </dl>
            </div>
            <div class="metric-card">
                <dl>
                    <dt>AFDB ID</dt>
                    <dd>ID of closest AlphaFold DB match. Click to view entry.</dd>
                </dl>
            </div>
        </div>

        <div class="toolbar">
            <div class="search-box">
                <input type="text" id="quickFilter" placeholder="Search all columns...">
            </div>
            <button class="btn btn-primary" onclick="exportToCsv()">Export CSV</button>
            <button class="btn btn-secondary" onclick="resetFilters()">Reset Filters</button>
        </div>

        <div id="metricsGrid" class="ag-theme-alpine"></div>
    </div>

    <!-- Sequence Modal -->
    <div id="sequenceModal" class="modal">
        <div class="modal-content">
            <div class="modal-header">
                <h3 id="modalTitle">Sequence</h3>
                <button class="close-btn" onclick="closeModal()">&times;</button>
            </div>
            <div id="modalSequence" class="sequence-display"></div>
            <button class="btn btn-secondary copy-btn" onclick="copySequence()">Copy to Clipboard</button>
        </div>
    </div>

    <!-- AG Grid JS -->
    <script src="https://cdn.jsdelivr.net/npm/ag-grid-community@32.3.3/dist/ag-grid-community.min.js"></script>

    <script>
        // Row data from Python
        const rowData = {rows_json};

        // Store current sequence for copying
        let currentSequence = '';

        // Store grid API reference
        let gridApi = null;

        // Custom cell renderer for sequence (collapsible)
        function sequenceCellRenderer(params) {{
            if (!params.value) return '-';

            const seq = params.value;
            const preview = seq.length > 15 ? seq.substring(0, 15) + '...' : seq;

            return `<span class="sequence-cell sequence-collapsed" onclick="showSequence('${{params.data.design_id}}', '${{seq}}')">${{preview}} <span class="expand-icon">&#x25BC;</span></span>`;
        }}

        // Custom cell renderer for AFDB link
        function afdbCellRenderer(params) {{
            if (!params.value) return '-';

            const url = params.data.afdb_url;
            if (url) {{
                return `<a href="${{url}}" target="_blank" class="afdb-link">${{params.value}}</a>`;
            }}
            return params.value;
        }}

        // Custom cell renderer for Foldseek target (clickable to copy)
        function foldseekTargetCellRenderer(params) {{
            if (!params.value) return '-';

            const target = params.value;
            const preview = target.length > 20 ? target.substring(0, 20) + '...' : target;

            return `<span class="sequence-cell sequence-collapsed" onclick="showFoldseekTarget('${{params.data.design_id}}', '${{target}}')">${{preview}} <span class="expand-icon">&#x25BC;</span></span>`;
        }}

        // Custom cell renderer for ipSAE with color coding
        function ipsaeCellRenderer(params) {{
            if (params.value === null || params.value === undefined) return '-';

            const val = params.value.toFixed(4);
            if (params.value < 0.5) {{
                return `<span class="good-value">${{val}}</span>`;
            }}
            return val;
        }}

        // Custom cell renderer for binding affinity with color coding
        function affinityCellRenderer(params) {{
            if (params.value === null || params.value === undefined) return '-';

            const val = params.value.toFixed(1);
            if (params.value < -10) {{
                return `<span class="good-value">${{val}}</span>`;
            }}
            return val;
        }}

        // Column definitions
        const columnDefs = [
            {{
                headerName: '#',
                valueGetter: 'node.rowIndex + 1',
                pinned: 'left',
                filter: false,
                sortable: false,
                minWidth: 60,
                maxWidth: 70,
                suppressSizeToFit: true
            }},
            {{
                headerName: 'Design ID',
                field: 'design_id',
                pinned: 'left',
                filter: 'agTextColumnFilter',
                sortable: true,
                minWidth: 140
            }},
            {{
                headerName: 'ipSAE',
                field: 'ipsae',
                filter: 'agNumberColumnFilter',
                sortable: true,
                cellRenderer: ipsaeCellRenderer,
                minWidth: 100
            }},
            {{
                headerName: 'Binding Affinity (kcal/mol)',
                field: 'binding_affinity',
                filter: 'agNumberColumnFilter',
                sortable: true,
                cellRenderer: affinityCellRenderer,
                minWidth: 160
            }},
            {{
                headerName: 'AA Length',
                field: 'aa_length',
                filter: 'agNumberColumnFilter',
                sortable: true,
                minWidth: 100,
                valueFormatter: params => params.value ? params.value : '-'
            }},
            {{
                headerName: 'Sequence',
                field: 'sequence',
                filter: 'agTextColumnFilter',
                sortable: false,
                cellRenderer: sequenceCellRenderer,
                minWidth: 150
            }},
            {{
                headerName: 'AFDB ID',
                field: 'afdb_id',
                filter: 'agTextColumnFilter',
                sortable: true,
                cellRenderer: afdbCellRenderer,
                minWidth: 120
            }},
            {{
                headerName: 'Foldseek Distance',
                field: 'foldseek_distance',
                filter: 'agNumberColumnFilter',
                sortable: true,
                minWidth: 140,
                valueFormatter: params => params.value !== null && params.value !== undefined ? params.value : '-'
            }},
            {{
                headerName: 'Foldseek Target',
                field: 'foldseek_target',
                filter: 'agTextColumnFilter',
                sortable: true,
                minWidth: 200,
                cellRenderer: foldseekTargetCellRenderer
            }}
        ];

        // Grid options
        const gridOptions = {{
            columnDefs: columnDefs,
            rowData: rowData,
            defaultColDef: {{
                resizable: true,
                floatingFilter: true
            }},
            animateRows: true,
            pagination: true,
            paginationPageSize: 25,
            paginationPageSizeSelector: [10, 25, 50, 100],
            domLayout: 'normal'
        }};

        // Initialize grid
        document.addEventListener('DOMContentLoaded', function() {{
            const gridDiv = document.querySelector('#metricsGrid');
            gridApi = agGrid.createGrid(gridDiv, gridOptions);

            // Quick filter
            document.getElementById('quickFilter').addEventListener('input', function(e) {{
                gridApi.setGridOption('quickFilterText', e.target.value);
            }});
        }});

        // Export to CSV
        function exportToCsv() {{
            gridApi.exportDataAsCsv({{
                fileName: 'protein_design_metrics.csv',
                columnKeys: ['design_id', 'ipsae', 'binding_affinity', 'aa_length', 'sequence', 'afdb_id', 'foldseek_distance', 'foldseek_target']
            }});
        }}

        // Reset filters
        function resetFilters() {{
            gridApi.setFilterModel(null);
            document.getElementById('quickFilter').value = '';
            gridApi.setGridOption('quickFilterText', '');
        }}

        // Show sequence modal
        function showSequence(designId, sequence) {{
            currentSequence = sequence;
            document.getElementById('modalTitle').textContent = `Sequence: ${{designId}} (${{sequence.length}} AA)`;
            document.getElementById('modalSequence').textContent = sequence;
            document.getElementById('sequenceModal').style.display = 'block';
        }}

        // Show Foldseek target modal
        function showFoldseekTarget(designId, target) {{
            currentSequence = target;
            document.getElementById('modalTitle').textContent = `Foldseek Target: ${{designId}}`;
            document.getElementById('modalSequence').textContent = target;
            document.getElementById('sequenceModal').style.display = 'block';
        }}

        // Close modal
        function closeModal() {{
            document.getElementById('sequenceModal').style.display = 'none';
        }}

        // Copy sequence to clipboard
        function copySequence() {{
            navigator.clipboard.writeText(currentSequence).then(function() {{
                const btn = document.querySelector('.copy-btn');
                const originalText = btn.textContent;
                btn.textContent = 'Copied!';
                setTimeout(() => btn.textContent = originalText, 2000);
            }});
        }}

        // Close modal when clicking outside
        window.onclick = function(event) {{
            const modal = document.getElementById('sequenceModal');
            if (event.target === modal) {{
                closeModal();
            }}
        }}

        // Close modal with Escape key
        document.addEventListener('keydown', function(e) {{
            if (e.key === 'Escape') {{
                closeModal();
            }}
        }});
    </script>
</body>
</html>
'''

    with open(output_file, 'w') as f:
        f.write(html)

    print(f"HTML report generated: {output_file}")


def generate_csv_report(metrics, output_file):
    """Generate a CSV summary for programmatic access."""
    rows = []
    for design_id, data in sorted(metrics.items()):
        rows.append({
            'design_id': design_id,
            'ipsae': data.get('ipsae', ''),
            'binding_affinity_kcal_mol': data.get('binding_affinity', ''),
            'aa_length': data.get('aa_length', ''),
            'sequence': data.get('sequence', ''),
            'afdb_id': data.get('afdb_id', ''),
            'afdb_url': data.get('afdb_url', ''),
            'foldseek_target': data.get('foldseek_target', ''),
            'foldseek_distance': data.get('foldseek_distance', '')
        })

    # Sort by ipSAE ascending (lower = better)
    rows.sort(key=lambda x: x.get('ipsae') or float('inf'))

    fieldnames = ['design_id', 'ipsae', 'binding_affinity_kcal_mol', 'aa_length',
                  'sequence', 'afdb_id', 'afdb_url', 'foldseek_target', 'foldseek_distance']

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"CSV report generated: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Consolidate protein design metrics and generate HTML report"
    )
    parser.add_argument(
        "--ipsae_dir",
        default="ipsae",
        help="Directory containing ipSAE score files"
    )
    parser.add_argument(
        "--prodigy_dir",
        default="prodigy",
        help="Directory containing Prodigy result files"
    )
    parser.add_argument(
        "--foldseek_dir",
        default="foldseek",
        help="Directory containing Foldseek summary files"
    )
    parser.add_argument(
        "--sequence_dir",
        default=None,
        help="Directory containing binder sequence FASTA files"
    )
    parser.add_argument(
        "--output_html",
        default="design_metrics_report.html",
        help="Output HTML report filename"
    )
    parser.add_argument(
        "--output_csv",
        default="design_metrics_summary.csv",
        help="Output CSV summary filename"
    )
    parser.add_argument(
        "--title",
        default="Protein Design Metrics Report",
        help="Report title"
    )
    parser.add_argument(
        "--ipsae_cutoffs",
        default="10_10",
        help="ipSAE cutoffs suffix (e.g., '10_10' for PAE=10, dist=10)"
    )

    args = parser.parse_args()

    # Collect metrics from staged directories
    metrics = collect_metrics_from_dirs(
        args.ipsae_dir,
        args.prodigy_dir,
        args.foldseek_dir,
        args.sequence_dir,
        args.ipsae_cutoffs
    )

    if not metrics:
        print("Warning: No metrics found!", file=sys.stderr)
        # Create empty reports
        generate_html_report({}, args.output_html, args.title)
        generate_csv_report({}, args.output_csv)
        return

    print(f"\nCollected metrics for {len(metrics)} designs")

    # Generate reports
    generate_html_report(metrics, args.output_html, args.title)
    generate_csv_report(metrics, args.output_csv)


if __name__ == "__main__":
    main()
