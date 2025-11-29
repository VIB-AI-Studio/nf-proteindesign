process CONSOLIDATE_METRICS {
    label 'process_low'

    // Publish reports to top-level output directory
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    container 'community.wave.seqera.io/library/numpy:2.3.5--f8d2712d76b3e3ce'

    input:
    path ipsae_files, stageAs: 'ipsae/*'       // Collection of ipSAE score files
    path prodigy_files, stageAs: 'prodigy/*'   // Collection of Prodigy result files
    path foldseek_files, stageAs: 'foldseek/*' // Collection of Foldseek summary files
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

    # Debug: List staged files in subdirectories
    echo "=== Staged ipSAE files ==="
    ls -la ipsae/ 2>/dev/null || echo "No ipsae directory"
    echo ""
    echo "=== Staged Prodigy files ==="
    ls -la prodigy/ 2>/dev/null || echo "No prodigy directory"
    echo ""
    echo "=== Staged Foldseek files ==="
    ls -la foldseek/ 2>/dev/null || echo "No foldseek directory"
    echo ""

    # Run consolidation script with staged subdirectories
    python ${consolidate_script} \\
        --ipsae_dir "ipsae" \\
        --prodigy_dir "prodigy" \\
        --foldseek_dir "foldseek" \\
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
