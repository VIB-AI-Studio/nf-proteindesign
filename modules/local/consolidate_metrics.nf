process CONSOLIDATE_METRICS {
    label 'process_low'

    // Publish reports to top-level output directory
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    container 'community.wave.seqera.io/library/numpy:2.3.5--f8d2712d76b3e3ce'

    input:
    path ipsae_files, stageAs: 'ipsae/*'       // Collection of ipSAE score files
    path prodigy_files, stageAs: 'prodigy/*'   // Collection of Prodigy result files
    path foldseek_files, stageAs: 'foldseek/*' // Collection of Foldseek summary files
    path sequence_files, stageAs: 'sequences/*' // Collection of binder sequence files
    path consolidate_script

    output:
    path "design_metrics_report.html", emit: report_html
    path "design_metrics_summary.csv", emit: summary_csv
    path "versions.yml", emit: versions

    script:
    def pae_cutoff = params.ipsae_pae_cutoff ?: 10
    def dist_cutoff = params.ipsae_dist_cutoff ?: 10

    """
    # Make script executable
    chmod +x ${consolidate_script}

    # Create directories if they don't exist (handles empty input lists)
    mkdir -p ipsae prodigy foldseek sequences

    # Debug: List staged files in subdirectories
    echo "=== Staged ipSAE files ==="
    ls -la ipsae/ 2>/dev/null || echo "No ipsae files"
    echo ""
    echo "=== Staged Prodigy files ==="
    ls -la prodigy/ 2>/dev/null || echo "No prodigy files"
    echo ""
    echo "=== Staged Foldseek files ==="
    ls -la foldseek/ 2>/dev/null || echo "No foldseek files"
    echo ""
    echo "=== Staged Sequence files ==="
    ls -la sequences/ 2>/dev/null || echo "No sequence files"
    echo ""

    # Build command with optional sequence directory (only if files exist)
    SEQ_FLAG=""
    if [ -n "\$(ls -A sequences/ 2>/dev/null)" ]; then
        SEQ_FLAG="--sequence_dir sequences"
    fi

    # Run consolidation script with staged subdirectories
    python ${consolidate_script} \\
        --ipsae_dir "ipsae" \\
        --prodigy_dir "prodigy" \\
        --foldseek_dir "foldseek" \\
        \${SEQ_FLAG} \\
        --output_html design_metrics_report.html \\
        --output_csv design_metrics_summary.csv \\
        --title "Protein Design Metrics Report" \\
        --ipsae_cutoffs "${pae_cutoff}_${dist_cutoff}"

    # Generate version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch design_metrics_report.html
    touch design_metrics_summary.csv
    touch versions.yml
    """
}
