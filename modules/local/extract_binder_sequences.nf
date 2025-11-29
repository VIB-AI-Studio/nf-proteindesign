/*
========================================================================================
    EXTRACT_BINDER_SEQUENCES: Extract binder sequences from Boltz2 structures
========================================================================================
    This process extracts the binder chain sequence from Boltz2 CIF structures
    for inclusion in the final metrics report.
----------------------------------------------------------------------------------------
*/

process EXTRACT_BINDER_SEQUENCES {
    tag "${meta.id}"
    label 'process_low'

    container 'community.wave.seqera.io/library/pip_biopython:79c9733809036fb1'

    input:
    tuple val(meta), path(structure)

    output:
    tuple val(meta), path("${meta.id}_binder_sequence.fa"), emit: sequence
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import sys
    from pathlib import Path

    try:
        from Bio import PDB
        from Bio.PDB import MMCIFParser, PDBParser
        from Bio.SeqUtils import seq1
        import Bio
    except ImportError:
        print("ERROR: Biopython not installed", file=sys.stderr)
        sys.exit(1)

    structure_file = "${structure}"
    design_id = "${meta.id}"

    # Parse structure
    if structure_file.endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure('structure', structure_file)

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
                        one_letter = seq1(resname)
                        residues.append(one_letter)
                    except (KeyError, ValueError):
                        residues.append('X')

            if residues:
                sequences[chain_id] = ''.join(residues)

    if not sequences:
        print("ERROR: No amino acid sequences found", file=sys.stderr)
        sys.exit(1)

    # Identify binder chain (always the shortest chain)
    # In Boltz2 output, the binder is the designed smaller protein
    # The target is typically the larger protein we're designing against
    binder_chain_id = min(sequences.items(), key=lambda x: len(x[1]))[0]
    binder_sequence = sequences[binder_chain_id]

    print(f"Chains found: {{{', '.join(f'{k}: {len(v)} AA' for k, v in sequences.items())}}}")
    print(f"Selected binder chain {binder_chain_id} (shortest) with {len(binder_sequence)} AA")

    # Write FASTA output
    with open(f"{design_id}_binder_sequence.fa", 'w') as f:
        f.write(f">{design_id}:B binder_sequence length={len(binder_sequence)}\\n")
        f.write(binder_sequence + "\\n")

    print(f"Extracted binder sequence from chain {binder_chain_id}: {len(binder_sequence)} AA")

    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"NFPROTEINDESIGN:PROTEIN_DESIGN:EXTRACT_BINDER_SEQUENCES":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
        f.write(f'    biopython: {Bio.__version__}\\n')
    """

    stub:
    """
    echo ">stub_sequence" > ${meta.id}_binder_sequence.fa
    echo "MOCKSEQUENCE" >> ${meta.id}_binder_sequence.fa
    touch versions.yml
    """
}
