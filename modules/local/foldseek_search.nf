process FOLDSEEK_SEARCH {
    tag "${meta.id}"
    label 'process_medium'

    // Publish results
    publishDir "${params.outdir}/${meta.parent_id ?: meta.id}/foldseek", mode: params.publish_dir_mode

    container 'ghcr.io/steineggerlab/foldseek:master-cuda12'

    // GPU acceleration - Foldseek supports GPU for faster searches (4-27x speedup)
    accelerator 1, type: 'nvidia-gpu'

    input:
    tuple val(meta), path(structure)
    path database

    output:
    tuple val(meta), path("${meta.id}_foldseek_results.tsv"), emit: results
    tuple val(meta), path("${meta.id}_foldseek_summary.tsv"), emit: summary
    path "versions.yml", emit: versions

    script:
    // Determine database path - can be a path or directory
    def db_path = database.name != 'NO_DATABASE' ? database : params.foldseek_database

    // Set search parameters
    def evalue = params.foldseek_evalue ?: 0.001
    def max_seqs = params.foldseek_max_seqs ?: 100
    def sensitivity = params.foldseek_sensitivity ?: 9.5
    def coverage = params.foldseek_coverage ?: 0.0
    def alignment_type = params.foldseek_alignment_type ?: 2

    // Validate database
    if (!db_path) {
        error "ERROR: No Foldseek database specified. Please set --foldseek_database parameter."
    }

    """
    /usr/local/bin/foldseek_avx2 easy-search \\
        ${structure} \\
        ${db_path}/afdb \\
        ${meta.id}_foldseek_results.tsv \\
        tmp_foldseek \\
        -e ${evalue} \\
        --max-seqs ${max_seqs} \\
        -s ${sensitivity} \\
        -c ${coverage} \\
        --alignment-type ${alignment_type} \\
        --threads ${task.cpus} \\
        --gpu 1 \\
        --prefilter-mode 1

    # Create summary (top hits)
    head -20 ${meta.id}_foldseek_results.tsv > ${meta.id}_foldseek_summary.tsv

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$(foldseek version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_foldseek_results.tsv
    touch ${meta.id}_foldseek_summary.tsv
    touch versions.yml
    """
}
