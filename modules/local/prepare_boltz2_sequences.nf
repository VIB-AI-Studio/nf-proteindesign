/*
========================================================================================
    PREPARE_BOLTZ2_SEQUENCES: Prepare sequences for Boltz2 refolding
========================================================================================
    This process combines two operations:
    1. Extract target sequence from original Boltzgen structure (for Boltz2 input)
    2. Split ProteinMPNN multi-sequence FASTA into individual files
    
    All sequences are included (original Boltzgen sequence + MPNN-designed sequences).
----------------------------------------------------------------------------------------
*/

process PREPARE_BOLTZ2_SEQUENCES {
    tag "${meta.id}"
    label 'process_low'
    
    container 'community.wave.seqera.io/library/pip_biopython:79c9733809036fb1'

    input:
    tuple val(meta), path(mpnn_fasta), path(boltzgen_structure)

    output:
    tuple val(meta), path("sequences/*.fa"), emit: sequences
    tuple val(meta), path("${meta.id}_target_sequence.txt"), emit: target_sequence
    tuple val(meta), path("${meta.id}_target_info.json"), emit: target_info
    path "versions.yml", emit: versions

    script:
    """
    prepare_boltz2_sequences.py \\
        ${mpnn_fasta} \\
        ${boltzgen_structure} \\
        ${meta.id}
    """

    stub:
    """
    mkdir -p sequences
    echo "MOCK_TARGET_SEQUENCE" > ${meta.id}_target_sequence.txt
    echo '{"design_id": "${meta.id}", "target_chain": "A", "target_length": 100}' > ${meta.id}_target_info.json
    touch sequences/placeholder_seq_0.fa
    touch versions.yml
    """
}
