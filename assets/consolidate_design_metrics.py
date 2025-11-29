#!/usr/bin/env python3
"""
Consolidate protein design metrics from ipSAE, Prodigy, and Foldseek outputs.

Generates an interactive HTML report with a searchable, filterable table.

This script reads from staged files in a flat directory structure.
"""

import argparse
import csv
import os
import sys
import re
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

                    # Calculate distance
                    distance = alnlen * (1.0 - fident)

                    if distance < best_distance:
                        best_distance = distance
                        best_target = target
                except (ValueError, IndexError):
                    continue

        if best_target:
            return best_target, round(best_distance, 1)

    except Exception as e:
        print(f"Warning: Could not parse Foldseek file {foldseek_file}: {e}", file=sys.stderr)

    return None, None


def extract_design_id_from_ipsae(filename):
    """
    Extract design ID from ipSAE filename.

    Example: 2vsm_r1_s0_seq_0_model_0_10_10.txt -> 2vsm_r1_s0
    """
    base = Path(filename).stem
    # Remove _XX_XX cutoff pattern at the end
    base = re.sub(r'_\d+_\d+$', '', base)
    # Remove _seq_X_model_Y pattern
    base = re.sub(r'_seq_\d+_model_\d+$', '', base)
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


def collect_metrics_from_dirs(ipsae_dir, prodigy_dir, foldseek_dir, ipsae_cutoffs="10_10"):
    """
    Collect metrics from staged files in separate directories.

    Args:
        ipsae_dir: Directory containing ipSAE score files
        prodigy_dir: Directory containing Prodigy result files
        foldseek_dir: Directory containing Foldseek summary files
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
                print(f"  Foldseek [{design_id}]: {target} (dist={distance})")
    else:
        print(f"Foldseek directory not found: {foldseek_dir}")

    return metrics


def generate_html_report(metrics, output_file, title="Protein Design Metrics Report"):
    """
    Generate an interactive HTML report with DataTables for searching/filtering.
    """
    # Convert metrics dict to sorted list
    rows = []
    for design_id, data in sorted(metrics.items()):
        rows.append({
            'design_id': design_id,
            'ipsae': data.get('ipsae'),
            'binding_affinity': data.get('binding_affinity'),
            'foldseek_target': data.get('foldseek_target'),
            'foldseek_distance': data.get('foldseek_distance')
        })

    # Sort by ipSAE descending (higher = better interface prediction)
    rows.sort(key=lambda x: x.get('ipsae') or 0, reverse=True)

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>

    <!-- DataTables CSS -->
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.7/css/jquery.dataTables.min.css">
    <link rel="stylesheet" href="https://cdn.datatables.net/buttons/2.4.2/css/buttons.dataTables.min.css">

    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }}

        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}

        h1 {{
            color: #333;
            margin-bottom: 10px;
        }}

        .summary {{
            color: #666;
            margin-bottom: 25px;
            padding: 15px;
            background: #f8f9fa;
            border-radius: 5px;
        }}

        .metric-info {{
            font-size: 0.9em;
            color: #666;
            margin-bottom: 20px;
        }}

        .metric-info dt {{
            font-weight: bold;
            color: #333;
        }}

        .metric-info dd {{
            margin-left: 20px;
            margin-bottom: 8px;
        }}

        table.dataTable {{
            width: 100% !important;
        }}

        table.dataTable thead th {{
            background-color: #4a90d9;
            color: white;
        }}

        table.dataTable tbody tr:hover {{
            background-color: #e8f4f8 !important;
        }}

        .good-value {{
            color: #28a745;
            font-weight: bold;
        }}

        .medium-value {{
            color: #ffc107;
        }}

        .highlight-row {{
            background-color: #d4edda !important;
        }}

        .dt-buttons {{
            margin-bottom: 15px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>{title}</h1>

        <div class="summary">
            <strong>Total designs evaluated:</strong> {len(rows)}
        </div>

        <div class="metric-info">
            <dl>
                <dt>ipSAE (Interface Predicted Structural Alignment Error)</dt>
                <dd>Lower is better. Measures interface quality; values &lt; 0.5 indicate good interfaces.</dd>

                <dt>Binding Affinity (kcal/mol)</dt>
                <dd>More negative is better. Predicted from PRODIGY; values &lt; -10 indicate strong binding.</dd>

                <dt>Foldseek Distance</dt>
                <dd>AA differences to closest structural match. Lower indicates similarity to known structures.</dd>
            </dl>
        </div>

        <table id="metricsTable" class="display" style="width:100%">
            <thead>
                <tr>
                    <th>Design ID</th>
                    <th>ipSAE</th>
                    <th>Binding Affinity (kcal/mol)</th>
                    <th>Foldseek Target</th>
                    <th>Foldseek Distance</th>
                </tr>
            </thead>
            <tbody>
'''

    for row in rows:
        # Format values
        ipsae = f"{row['ipsae']:.4f}" if row['ipsae'] is not None else '-'
        affinity = f"{row['binding_affinity']:.1f}" if row['binding_affinity'] is not None else '-'
        fs_target = row['foldseek_target'] or '-'
        fs_dist = f"{row['foldseek_distance']:.1f}" if row['foldseek_distance'] is not None else '-'

        # Add CSS classes for good values
        ipsae_class = 'good-value' if row['ipsae'] is not None and row['ipsae'] < 0.5 else ''
        affinity_class = 'good-value' if row['binding_affinity'] is not None and row['binding_affinity'] < -10 else ''

        html += f'''                <tr>
                    <td><strong>{row['design_id']}</strong></td>
                    <td class="{ipsae_class}">{ipsae}</td>
                    <td class="{affinity_class}">{affinity}</td>
                    <td>{fs_target}</td>
                    <td>{fs_dist}</td>
                </tr>
'''

    html += '''            </tbody>
        </table>
    </div>

    <!-- jQuery -->
    <script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>

    <!-- DataTables JS -->
    <script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/2.4.2/js/dataTables.buttons.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/2.4.2/js/buttons.html5.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.10.1/jszip.min.js"></script>

    <script>
        $(document).ready(function() {
            $('#metricsTable').DataTable({
                pageLength: 25,
                order: [[1, 'desc']],  // Sort by ipSAE descending
                dom: 'Bfrtip',
                buttons: [
                    {
                        extend: 'csv',
                        text: 'Export CSV',
                        filename: 'protein_design_metrics'
                    },
                    {
                        extend: 'copy',
                        text: 'Copy to Clipboard'
                    }
                ],
                columnDefs: [
                    {
                        targets: [1, 2, 4],  // Numeric columns
                        type: 'num'
                    }
                ],
                language: {
                    search: "Filter designs:",
                    lengthMenu: "Show _MENU_ designs per page"
                }
            });
        });
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
            'foldseek_target': data.get('foldseek_target', ''),
            'foldseek_distance': data.get('foldseek_distance', '')
        })

    # Sort by ipSAE descending
    rows.sort(key=lambda x: x.get('ipsae') or 0, reverse=True)

    fieldnames = ['design_id', 'ipsae', 'binding_affinity_kcal_mol',
                  'foldseek_target', 'foldseek_distance']

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
