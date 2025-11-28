/*
========================================================================================
    MMSEQS2_MSA: Generate multiple sequence alignments using MMSeqs2 GPU
========================================================================================
    This process generates MSAs for protein sequences using MMSeqs2 with GPU acceleration.
    It creates MSAs in A3M format that can be used by Boltz-2 for improved predictions.
    
    Key features:
    - GPU-accelerated alignment using MMSeqs2
    - Deduplication to run MSA only once per unique sequence
    - Compatible with Boltz-2 MSA input format
----------------------------------------------------------------------------------------
*/

process MMSEQS2_MSA {
    tag "${meta.id}"
    label 'process_high_gpu'
    
    // Publish MSA results
    publishDir "${params.outdir}/${meta.parent_id}/msa", mode: params.publish_dir_mode

    // MMSeqs2 container with GPU support
    container 'quay.io/biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_2'
    
    // GPU acceleration for fast MSA generation
    accelerator 1, type: 'nvidia-gpu'

    input:
    tuple val(meta), path(sequence_fasta)

    output:
    tuple val(meta), path("${meta.id}_msa.a3m"), emit: msa
    tuple val(meta), path("${meta.id}_msa_stats.txt"), emit: stats
    path "versions.yml", emit: versions

    script:
    def database = params.mmseqs2_database ?: "\$MMSEQS2_DB"
    def evalue = params.mmseqs2_evalue ?: 1e-3
    def iterations = params.mmseqs2_iterations ?: 3
    def sensitivity = params.mmseqs2_sensitivity ?: 7.5
    def max_seqs = params.mmseqs2_max_seqs ?: 1000
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "============================================"
    echo "MMSeqs2 GPU-Accelerated MSA Generation"
    echo "============================================"
    
    # Check for GPU
    if command -v nvidia-smi &> /dev/null && nvidia-smi &> /dev/null; then
        echo "✓ GPU detected"
        nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
        USE_GPU="--gpu 1"
    else
        echo "⚠  No GPU detected - falling back to CPU (slower)"
        USE_GPU=""
    fi
    
    # Check if database is set
    if [ -z "${database}" ] || [ ! -d "${database}" ]; then
        echo "ERROR: MMSeqs2 database not specified or does not exist"
        echo "Please set params.mmseqs2_database to a valid database path"
        echo "Example databases: UniRef30, UniRef90, ColabFoldDB"
        exit 1
    fi
    
    echo ""
    echo "Configuration:"
    echo "  Database: ${database}"
    echo "  E-value threshold: ${evalue}"
    echo "  Iterations: ${iterations}"
    echo "  Sensitivity: ${sensitivity}"
    echo "  Max sequences: ${max_seqs}"
    echo "  GPU acceleration: \${USE_GPU:-disabled}"
    
    # Create temporary directory
    mkdir -p tmp_msa
    
    # Create MMSeqs2 query database from input FASTA
    echo ""
    echo "Creating query database..."
    mmseqs createdb ${sequence_fasta} queryDB
    
    # Prepare target database for GPU (if using GPU)
    if [ -n "\${USE_GPU}" ]; then
        echo ""
        echo "Preparing database for GPU acceleration..."
        
        # Check if database is already in padded format
        if [ ! -f "${database}.dbtype" ]; then
            echo "Creating padded sequence database for GPU..."
            mmseqs makepaddedseqdb ${database} ${database}_gpu tmp_msa
            TARGET_DB="${database}_gpu"
        else
            # Check if it's already a padded database
            DB_TYPE=\$(cat ${database}.dbtype 2>/dev/null || echo "0")
            if [ "\${DB_TYPE}" = "15" ]; then
                echo "Database is already in GPU-compatible format"
                TARGET_DB="${database}"
            else
                echo "Creating padded sequence database for GPU..."
                mmseqs makepaddedseqdb ${database} targetDB_gpu tmp_msa
                TARGET_DB="targetDB_gpu"
            fi
        fi
    else
        TARGET_DB="${database}"
    fi
    
    # Run MMSeqs2 search with GPU acceleration
    echo ""
    echo "Running MMSeqs2 easy-search..."
    mmseqs easy-search \\
        ${sequence_fasta} \\
        \${TARGET_DB} \\
        search_results.m8 \\
        tmp_msa \\
        -e ${evalue} \\
        --num-iterations ${iterations} \\
        -s ${sensitivity} \\
        --max-seqs ${max_seqs} \\
        --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits" \\
        \${USE_GPU}
    
    # Convert search results to A3M format (MSA format for Boltz-2)
    echo ""
    echo "Converting results to A3M format..."
    
    # Create result database for MSA conversion
    mmseqs createdb ${sequence_fasta} queryDB_msa
    mmseqs search queryDB_msa \${TARGET_DB} resultDB tmp_msa \\
        -e ${evalue} \\
        --num-iterations ${iterations} \\
        -s ${sensitivity} \\
        --max-seqs ${max_seqs} \\
        \${USE_GPU}
    
    # Convert to MSA format (A3M)
    mmseqs result2msa queryDB_msa \${TARGET_DB} resultDB msaDB \\
        --msa-format-mode 6
    
    # Export to A3M file
    mmseqs msa2profile msaDB msaDB_prof --match-mode 1
    mmseqs result2profile msaDB_prof msaDB_prof
    mmseqs convertalis queryDB_msa \${TARGET_DB} resultDB ${meta.id}_msa.a3m \\
        --format-mode 3
    
    # If A3M creation failed, try alternative method
    if [ ! -s "${meta.id}_msa.a3m" ]; then
        echo "Primary A3M conversion failed, trying alternative method..."
        
        # Create A3M from search results
        mmseqs convertalis queryDB_msa \${TARGET_DB} resultDB ${meta.id}_msa.a3m \\
            --format-mode 3 \\
            --format-output query,target,qaln,taln
    fi
    
    # If still failed, create minimal A3M with query sequence only
    if [ ! -s "${meta.id}_msa.a3m" ]; then
        echo "⚠  Warning: Could not create MSA, using query sequence only"
        cat ${sequence_fasta} > ${meta.id}_msa.a3m
    fi
    
    # Generate statistics
    echo ""
    echo "Generating MSA statistics..."
    
    # Count sequences in MSA
    MSA_COUNT=\$(grep -c "^>" ${meta.id}_msa.a3m || echo "0")
    
    # Calculate MSA depth and coverage
    python3 <<STATS
import os

# Read A3M file
with open('${meta.id}_msa.a3m', 'r') as f:
    lines = f.readlines()

# Parse sequences
sequences = []
current_seq = []
for line in lines:
    if line.startswith('>'):
        if current_seq:
            sequences.append(''.join(current_seq))
            current_seq = []
    else:
        current_seq.append(line.strip())
if current_seq:
    sequences.append(''.join(current_seq))

# Calculate statistics
num_seqs = len(sequences)
if num_seqs > 0:
    query_len = len(sequences[0].replace('-', ''))
    avg_len = sum(len(seq.replace('-', '')) for seq in sequences) / num_seqs
    coverage = (num_seqs - 1) / max(query_len, 1)  # MSA depth relative to query length
else:
    query_len = 0
    avg_len = 0
    coverage = 0

# Write statistics
with open('${meta.id}_msa_stats.txt', 'w') as f:
    f.write(f"MMSeqs2 MSA Statistics\\n")
    f.write(f"=====================\\n\\n")
    f.write(f"Sample ID: ${meta.id}\\n")
    f.write(f"Parent ID: ${meta.parent_id}\\n\\n")
    f.write(f"Query sequence length: {query_len}\\n")
    f.write(f"Number of sequences in MSA: {num_seqs}\\n")
    f.write(f"Average sequence length: {avg_len:.1f}\\n")
    f.write(f"MSA depth (sequences per residue): {coverage:.2f}\\n\\n")
    f.write(f"Parameters:\\n")
    f.write(f"  Database: ${database}\\n")
    f.write(f"  E-value: ${evalue}\\n")
    f.write(f"  Iterations: ${iterations}\\n")
    f.write(f"  Sensitivity: ${sensitivity}\\n")
    f.write(f"  Max sequences: ${max_seqs}\\n")
    f.write(f"  GPU: {'Yes' if '${USE_GPU:-}' else 'No'}\\n")

print(f"\\nMSA Statistics:")
print(f"  Total sequences: {num_seqs}")
print(f"  Query length: {query_len}")
print(f"  MSA depth: {coverage:.2f}")
STATS
    
    echo ""
    echo "============================================"
    echo "MMSeqs2 MSA Generation Complete"
    echo "============================================"
    echo "MSA file: ${meta.id}_msa.a3m"
    echo "Statistics: ${meta.id}_msa_stats.txt"
    echo "Sequences in MSA: \${MSA_COUNT}"
    echo "============================================"
    
    # Clean up temporary files
    rm -rf tmp_msa queryDB* resultDB* msaDB* targetDB_gpu* 2>/dev/null || true
    
    # Generate version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs2: \$(mmseqs version 2>&1 | head -n1 | awk '{print \$2}')
        python: \$(python3 --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_msa.a3m
    echo "Stub MSA statistics" > ${meta.id}_msa_stats.txt
    touch versions.yml
    """
}
