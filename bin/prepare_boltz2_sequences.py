#!/usr/bin/env python3
"""
Prepare sequences for Boltz2 refolding.

This script combines two operations:
1. Extract target sequence from original Boltzgen structure (for Boltz2 input)
2. Split ProteinMPNN multi-sequence FASTA into individual files

All sequences are included (original Boltzgen sequence + MPNN-designed sequences).
"""

import os
import sys
import json
from pathlib import Path
import argparse
import Bio
from Bio import PDB
from Bio.PDB import PDBIO, MMCIFParser, PDBParser
from Bio.SeqUtils import seq1


def extract_target_sequence(structure_file, design_id):
    """Extract target sequence from Boltzgen structure."""
    print("=" * 80)
    print("PART 1: Extracting target sequence from Boltzgen structure")
    print("=" * 80)
    
    structure_path = Path(structure_file)
    print(f"Structure file: {structure_path}")
    
    try:
        if structure_path.suffix.lower() == '.cif':
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)
        
        structure = parser.get_structure('structure', str(structure_path))
        
        # Extract sequences from all chains
        sequences = {}
        for model in structure:
            for chain in model:
                chain_id = chain.id
                residues = []
                for residue in chain:
                    if PDB.is_aa(residue):
                        resname = residue.get_resname()
                        try:
                            # seq1 converts 3-letter to 1-letter code
                            one_letter = seq1(resname)
                            residues.append(one_letter)
                        except (KeyError, ValueError):
                            # Handle non-standard amino acids
                            residues.append('X')
                
                if residues:
                    sequences[chain_id] = ''.join(residues)
        
        if not sequences:
            print("ERROR: No amino acid sequences found in structure", file=sys.stderr)
            sys.exit(1)
        
        # Identify target chain (longest chain)
        target_chain_id = max(sequences.items(), key=lambda x: len(x[1]))[0]
        target_sequence = sequences[target_chain_id]
        
        # Write target sequence to file
        target_seq_file = f"{design_id}_target_sequence.txt"
        with open(target_seq_file, 'w') as f:
            f.write(target_sequence + "\n")
        
        print(f"✓ Target chain {target_chain_id} extracted ({len(target_sequence)} residues)")
        
        # Create info JSON
        info = {
            "design_id": design_id,
            "source_structure": structure_path.name,
            "target_chain": target_chain_id,
            "target_length": len(target_sequence),
            "all_chains": {cid: len(seq) for cid, seq in sequences.items()}
        }
        
        target_info_file = f"{design_id}_target_info.json"
        with open(target_info_file, 'w') as f:
            json.dump(info, f, indent=2)
        
        return target_seq_file, target_info_file
        
    except Exception as e:
        print(f"ERROR: Failed to extract target sequence: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


def split_fasta(input_file, output_dir="sequences"):
    """Split multi-sequence FASTA into individual files."""
    print()
    print("=" * 80)
    print("PART 2: Splitting ProteinMPNN sequences")
    print("=" * 80)
    
    print(f"FASTA file: {input_file}")
    
    sequences = []
    current_seq = []
    current_header = None
    
    # Read FASTA
    try:
        with open(input_file, 'r') as f:
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
            
            # Add last sequence
            if current_header and current_seq:
                sequences.append((current_header, ''.join(current_seq)))
    except Exception as e:
        print(f"ERROR: Failed to read FASTA file {input_file}: {e}", file=sys.stderr)
        sys.exit(1)
        
    print(f"Found {len(sequences)} sequences in {input_file}")
    
    # Include ALL sequences (including the original first sequence)
    sequences_to_process = sequences
    
    if not sequences_to_process:
        print("ERROR: No sequences found in FASTA file", file=sys.stderr)
        sys.exit(1)
        
    print(f"Splitting {len(sequences_to_process)} sequences (including original)")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Write each sequence to a separate file
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    
    output_files = []
    for idx, (header, seq) in enumerate(sequences_to_process):
        # Use 0-based indexing: seq_0 is original, seq_1+ are MPNN designs
        seq_num = idx
        output_file = os.path.join(output_dir, f"{base_name}_seq_{seq_num}.fa")
        
        with open(output_file, 'w') as out:
            out.write(f"{header}\n{seq}\n")
        
        output_files.append(output_file)
        seq_type = "original" if idx == 0 else f"MPNN design {idx}"
        print(f"✓ Created {output_file} ({seq_type})")
    
    return output_files


def main():
    parser = argparse.ArgumentParser(description="Prepare sequences for Boltz2 refolding")
    parser.add_argument("mpnn_fasta", help="ProteinMPNN multi-sequence FASTA file")
    parser.add_argument("boltzgen_structure", help="Boltzgen structure file (CIF or PDB)")
    parser.add_argument("design_id", help="Design ID for output files")
    
    args = parser.parse_args()
    
    # Extract target sequence
    target_seq_file, target_info_file = extract_target_sequence(
        args.boltzgen_structure,
        args.design_id
    )
    
    # Split MPNN sequences
    sequence_files = split_fasta(args.mpnn_fasta)
    
    # Generate version information
    with open("versions.yml", "w") as f:
        f.write('"NFPROTEINDESIGN:PROTEIN_DESIGN:PREPARE_BOLTZ2_SEQUENCES":\n')
        f.write(f'    python: {sys.version.split()[0]}\n')
        f.write(f'    biopython: {Bio.__version__}\n')
    
    print()
    print("=" * 80)
    print("✓ Sequence preparation complete")
    print("=" * 80)


if __name__ == "__main__":
    main()