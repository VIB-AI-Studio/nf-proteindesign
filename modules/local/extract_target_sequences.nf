process EXTRACT_TARGET_SEQUENCES {
    tag "${meta.id}"
    label 'process_low'
    
    container 'biopython/biopython:latest'

    input:
    tuple val(meta), path(original_structures)

    output:
    tuple val(meta), path("${meta.id}_target_sequences.txt"), emit: target_sequences
    tuple val(meta), path("${meta.id}_target_info.json"), emit: target_info
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import os
    import sys
    import json
    from pathlib import Path
    import Bio
    from Bio import PDB
    from Bio.PDB import PDBIO, MMCIFParser, PDBParser
    
    # Find structure files
    structures_input = Path("${original_structures}")
    
    if structures_input.is_dir():
        structure_files = list(structures_input.rglob("*.cif")) + list(structures_input.rglob("*.pdb"))
    elif structures_input.is_file():
        structure_files = [structures_input]
    else:
        structure_files = [Path(f) for f in "${original_structures}".split() if Path(f).exists()]
    
    print("Found " + str(len(structure_files)) + " structure files")
    
    if len(structure_files) == 0:
        print("ERROR: No structure files found", file=sys.stderr)
        sys.exit(1)
    
    # Use the first structure file to extract target sequence
    first_structure = structure_files[0]
    print("Extracting target sequence from: " + str(first_structure))
    
    # Parse structure using BioPython
    try:
        if first_structure.suffix.lower() == '.cif':
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)
        
        structure = parser.get_structure('structure', str(first_structure))
        
        # Extract sequences from all chains
        sequences = {}
        for model in structure:
            for chain in model:
                chain_id = chain.id
                residues = []
                for residue in chain:
                    if PDB.is_aa(residue):
                        # Get 3-letter code and convert to 1-letter
                        resname = residue.get_resname()
                        try:
                            one_letter = PDB.Polypeptide.three_to_one(resname)
                            residues.append(one_letter)
                        except KeyError:
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
        with open("${meta.id}_target_sequences.txt", 'w') as f:
            f.write(target_sequence + "\\n")
        
        print("Target chain " + target_chain_id + " extracted (" + str(len(target_sequence)) + " residues)")
        
        # Create info JSON
        info = {
            "design_id": "${meta.id}",
            "source_structure": str(first_structure.name),
            "target_chain": target_chain_id,
            "target_length": len(target_sequence),
            "num_structures": len(structure_files),
            "all_chains": {cid: len(seq) for cid, seq in sequences.items()}
        }
        
        with open("${meta.id}_target_info.json", 'w') as f:
            json.dump(info, f, indent=2)
        
        print("Target info saved to ${meta.id}_target_info.json")
        
    except Exception as e:
        print("ERROR: Failed to extract sequence: " + str(e), file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Generate version information
    with open("versions.yml", "w") as f:
        f.write("\\"${task.process}\\":\\n")
        f.write("    python: " + sys.version.split()[0] + "\\n")
        f.write("    biopython: " + Bio.__version__ + "\\n")
    """

    stub:
    """
    echo "MOCK_TARGET_SEQUENCE" > ${meta.id}_target_sequences.txt
    echo '{"design_id": "${meta.id}", "target_chain": "A", "target_length": 100}' > ${meta.id}_target_info.json
    touch versions.yml
    """
}
