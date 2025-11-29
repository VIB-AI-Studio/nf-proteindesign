/*
========================================================================================
    PROTEIN_DESIGN: Workflow for protein design using YAML specifications
========================================================================================
    This workflow uses pre-made design YAML files for protein design with Boltzgen
    and optional analysis modules.
----------------------------------------------------------------------------------------
*/
include { BOLTZGEN_RUN } from '../modules/local/boltzgen_run'
include { CONVERT_CIF_TO_PDB } from '../modules/local/convert_cif_to_pdb'
include { PROTEINMPNN_OPTIMIZE } from '../modules/local/proteinmpnn_optimize'
include { PREPARE_BOLTZ2_SEQUENCES } from '../modules/local/prepare_boltz2_sequences'
include { BOLTZ2_REFOLD } from '../modules/local/boltz2_refold'
include { IPSAE_CALCULATE } from '../modules/local/ipsae_calculate'
include { PRODIGY_PREDICT } from '../modules/local/prodigy_predict'
include { FOLDSEEK_SEARCH } from '../modules/local/foldseek_search'
include { CONSOLIDATE_METRICS } from '../modules/local/consolidate_metrics'

workflow PROTEIN_DESIGN {
    
    take:
    ch_input    // channel: [meta, design_yaml, structure_files, target_msa]
    ch_cache    // channel: path to cache directory or EMPTY_CACHE placeholder

    main:
    
    // ========================================================================
    // Run Boltzgen on design YAMLs
    // ========================================================================
    
    // Prepare Boltzgen input by removing target_msa (not needed for Boltzgen)
    ch_boltzgen_input = ch_input
        .map { meta, design_yaml, structure_files, target_msa ->
            [meta, design_yaml, structure_files]
        }
    
    // Run Boltzgen for each design in parallel
    BOLTZGEN_RUN(ch_boltzgen_input, ch_cache)
    
    // ========================================================================
    // ProteinMPNN: Optimize sequences for designed structures
    // ========================================================================
    if (params.run_proteinmpnn) {
        // Step 1: Convert CIF structures to PDB format (ProteinMPNN requires PDB)
        // Use budget_design_cifs which contains ONLY the budget designs (e.g., 2 structures if budget=2)
        // NOT all designs from results directory
        CONVERT_CIF_TO_PDB(BOLTZGEN_RUN.out.budget_design_cifs)
        
        // Step 2: Parallelize ProteinMPNN - run separately for each budget design
        // Use flatMap to create individual tasks per PDB file (one per budget iteration)
        ch_pdb_per_design = CONVERT_CIF_TO_PDB.out.pdb_files_all
            .flatMap { meta, pdb_files ->
                // Convert to list if single file and create defensive copy
                def pdb_list = pdb_files instanceof List ? new ArrayList(pdb_files) : [pdb_files]

                // Create a separate channel entry for each PDB file
                pdb_list.collect { pdb_file ->
                    // Extract rank number from filename (e.g., "rank1_2VSM_protein_design_1" -> "1")
                    def rank_num = pdb_file.baseName.replaceAll(/^rank(\d+)_.*/, '$1')

                    // Simplified naming: {sample}_r{rank}
                    def design_meta = [
                        id: "${meta.id}_r${rank_num}",
                        parent_id: meta.id,
                        rank_num: rank_num,
                        design_name: pdb_file.baseName  // Keep original for reference
                    ]

                    [design_meta, pdb_file]
                }
            }
        
        // Run ProteinMPNN on each design individually (parallel execution per budget design)
        PROTEINMPNN_OPTIMIZE(ch_pdb_per_design)
        
        // Use ProteinMPNN optimized structures for downstream analyses
        ch_final_designs_for_analysis = PROTEINMPNN_OPTIMIZE.out.optimized_designs
        
        // ====================================================================
        // Step 3: Prepare sequences for Boltz-2 refolding
        // ====================================================================
        // This combines:
        // 1. Extracting target sequence from Boltzgen structure
        // 2. Splitting ProteinMPNN FASTA into individual sequence files
        // ====================================================================
        if (params.run_boltz2_refold) {
            // Prepare input: combine MPNN FASTA files with Boltzgen structures
            // Use FIRST budget design CIF for target extraction (target is same across all designs)
            ch_prepare_input = PROTEINMPNN_OPTIMIZE.out.sequences
                .flatMap { meta, fasta_files ->
                    def fasta_list = fasta_files instanceof List ? new ArrayList(fasta_files) : [fasta_files]
                    fasta_list.collect { fasta_file ->
                        [meta, fasta_file]
                    }
                }
                .map { meta, fasta ->
                    [meta.parent_id, meta, fasta]
                }
                .combine(
                    BOLTZGEN_RUN.out.budget_design_cifs.map { meta, cif_files ->
                        def cif_list = cif_files instanceof List ? cif_files : [cif_files]
                        def first_cif = cif_list.sort()[0]
                        [meta.id, first_cif]
                    },
                    by: 0
                )
                .map { parent_id, meta, fasta, cif ->
                    [meta, fasta, cif]
                }
            
            // Run combined preparation process
            PREPARE_BOLTZ2_SEQUENCES(ch_prepare_input)
            
            // ================================================================
            // Prepare Target MSA from Samplesheet
            // ================================================================
            ch_target_msa = ch_input
                .map { meta, design_yaml, structure_files, target_msa ->
                    def msa_file = target_msa ?: file('NO_MSA')
                    [meta.id, msa_file]
                }
            
            // ================================================================
            // Create channel for Boltz-2 refolding
            // ================================================================
            ch_boltz2_input = PREPARE_BOLTZ2_SEQUENCES.out.sequences
                .flatMap { meta, fasta_files ->
                    def fasta_list = fasta_files instanceof List ? new ArrayList(fasta_files) : [fasta_files]
                    fasta_list.collect { fasta_file ->
                        def seq_num = fasta_file.baseName.replaceAll(/.*_seq_(\d+)$/, '$1')

                        // Simplified naming: {sample}_r{rank}_s{seq}
                        def seq_meta = [
                            id: "${meta.id}_s${seq_num}",
                            parent_id: meta.parent_id,
                            rank_num: meta.rank_num,
                            seq_num: seq_num,
                            mpnn_parent_id: meta.id,
                            sequence_name: fasta_file.baseName
                        ]
                        [seq_meta, fasta_file]
                    }
                }
                .map { meta, fasta -> 
                    [meta.mpnn_parent_id, meta, fasta]
                }
                .combine(
                    PREPARE_BOLTZ2_SEQUENCES.out.target_sequence.map { meta, seq -> 
                        [meta.id, seq]
                    },
                    by: 0
                )
                .map { mpnn_parent_id, meta, fasta, target_seq ->
                    [meta.parent_id, meta, fasta, target_seq]
                }
                .combine(ch_target_msa, by: 0)
                .map { parent_id, meta, fasta, target_seq, target_msa ->
                    [meta, fasta, target_seq, target_msa]
                }
            
            // Run Boltz-2 structure prediction with target MSA
            // NOTE: Boltz-2 will automatically add missing MSA info to binder
            // NOTE: Boltz-2 outputs NPZ files natively - no conversion needed!
            BOLTZ2_REFOLD(ch_boltz2_input)
        }
    } else {
        // Use Boltzgen outputs directly if ProteinMPNN is disabled
        ch_final_designs_for_analysis = BOLTZGEN_RUN.out.results
    }
    
    // ========================================================================
    // OPTIONAL: IPSAE scoring if enabled
    // ========================================================================
    // NOTE: IPSAE requires NPZ confidence files. We now support both:
    //   1. Boltzgen budget designs (native NPZ output)
    //   2. Boltz-2 refolded structures (native NPZ output - no conversion needed!)
    if (params.run_ipsae) {
        // Prepare IPSAE script as a value channel (reusable across all tasks)
        ch_ipsae_script = Channel.fromPath("${projectDir}/assets/ipsae.py", checkIfExists: true).first()
        
        // ====================================================================
        // Process Boltz-2 refolded structures
        // ====================================================================
        if (params.run_proteinmpnn && params.run_boltz2_refold) {
            // Get NPZ and CIF pairs directly from Boltz-2 (native NPZ output!)
            // Only use the best model (model_0) from each Boltz2 prediction
            ch_ipsae_input = BOLTZ2_REFOLD.out.structures
                .join(BOLTZ2_REFOLD.out.pae_npz, by: 0)
                .flatMap { meta, cif_files, npz_files ->
                    // Convert to lists if single files
                    def cif_list = cif_files instanceof List ? cif_files : [cif_files]
                    def npz_list = npz_files instanceof List ? npz_files : [npz_files]

                    // Filter to only model_0 (best model)
                    cif_list = cif_list.findAll { it.name.endsWith('model_0.cif') }
                    npz_list = npz_list.findAll { it.name.contains('model_0') }

                    // Create a map of basenames for matching
                    def npz_map = [:]
                    npz_list.each { npz_file ->
                        // Extract base name (without pae_ prefix)
                        def base_name = npz_file.baseName.replaceAll(/^pae_/, '')
                        npz_map[base_name] = npz_file
                    }

                    // Match CIF files with their NPZ files
                    cif_list.collect { cif_file ->
                        def base_name = cif_file.baseName
                        def npz_file = npz_map[base_name]

                        if (npz_file) {
                            // Use simplified naming from Boltz2 meta directly
                            def ipsae_meta = [
                                id: meta.id,
                                parent_id: meta.parent_id,
                                rank_num: meta.rank_num,
                                seq_num: meta.seq_num,
                                source: "boltz2"
                            ]

                            [ipsae_meta, npz_file, cif_file]
                        } else {
                            log.warn "⚠️  No matching NPZ file found for ${cif_file.name}"
                            null
                        }
                    }.findAll { it != null }
                }

            // Run IPSAE calculation
            IPSAE_CALCULATE(ch_ipsae_input, ch_ipsae_script)
        } else {
            log.warn "⚠️  IPSAE requested but ProteinMPNN/Boltz2 not enabled. Skipping IPSAE."
        }
    }
    
    // ========================================================================
    // OPTIONAL: PRODIGY binding affinity prediction if enabled
    // ========================================================================
    if (params.run_prodigy) {
        // Prepare PRODIGY parser script as a value channel (reusable across all tasks)
        ch_prodigy_script = Channel.fromPath("${projectDir}/assets/parse_prodigy_output.py", checkIfExists: true).first()
        
        if (params.run_proteinmpnn && params.run_boltz2_refold) {
            // Only use the best model (model_0) from each Boltz2 prediction
            ch_prodigy_input = BOLTZ2_REFOLD.out.structures
                .flatMap { meta, cif_files ->
                    // Convert to list if single file and create defensive copy
                    def cif_list = cif_files instanceof List ? new ArrayList(cif_files) : [cif_files]

                    // Filter to only model_0 (best model)
                    cif_list = cif_list.findAll { it.name.endsWith('model_0.cif') }

                    // Create a separate entry for each CIF file
                    cif_list.collect { cif_file ->
                        // Use simplified naming from Boltz2 meta directly
                        def design_meta = [
                            id: meta.id,
                            parent_id: meta.parent_id,
                            rank_num: meta.rank_num,
                            seq_num: meta.seq_num,
                            source: "boltz2"
                        ]

                        [design_meta, cif_file]
                    }
                }

            // Run PRODIGY binding affinity prediction
            PRODIGY_PREDICT(ch_prodigy_input, ch_prodigy_script)
        } else {
            log.warn "⚠️  Prodigy requested but ProteinMPNN/Boltz2 not enabled. Skipping Prodigy."
        }
    }
    
    // ========================================================================
    // OPTIONAL: Foldseek structural similarity search if enabled
    // ========================================================================
    // Search for structural homologs of both Boltzgen and Protenix structures
    // in the AlphaFold database (or other specified database)
    if (params.run_foldseek) {
        // Prepare database channel
        if (params.foldseek_database) {
            ch_foldseek_database = Channel.fromPath(params.foldseek_database, checkIfExists: true).first()
        } else {
            log.warn "⚠️  Foldseek is enabled but no database specified. Please set --foldseek_database parameter."
            ch_foldseek_database = Channel.value(file('NO_DATABASE'))
        }
        
        // ====================================================================
        // Process Boltz-2 refolded structures
        // ====================================================================
        if (params.run_proteinmpnn && params.run_boltz2_refold) {
            // Only use the best model (model_0) from each Boltz2 prediction
            ch_foldseek_input = BOLTZ2_REFOLD.out.structures
                .flatMap { meta, cif_files ->
                    // Convert to list if single file and create defensive copy
                    def cif_list = cif_files instanceof List ? new ArrayList(cif_files) : [cif_files]

                    // Filter to only model_0 (best model)
                    cif_list = cif_list.findAll { it.name.endsWith('model_0.cif') }

                    // Create a separate entry for each CIF file
                    cif_list.collect { cif_file ->
                        // Use simplified naming from Boltz2 meta directly
                        def design_meta = [
                            id: meta.id,
                            parent_id: meta.parent_id,
                            rank_num: meta.rank_num,
                            seq_num: meta.seq_num,
                            source: "boltz2"
                        ]

                        [design_meta, cif_file]
                    }
                }

            // Run Foldseek structural search
            FOLDSEEK_SEARCH(ch_foldseek_input, ch_foldseek_database)
        } else {
            log.warn "⚠️  Foldseek requested but ProteinMPNN/Boltz2 not enabled. Skipping Foldseek."
        }
    }
    
    // ========================================================================
    // CONSOLIDATION: Generate comprehensive metrics report
    // ========================================================================
    if (params.run_consolidation) {
        // Prepare consolidation script as a value channel (reusable)
        ch_consolidate_script = Channel.fromPath("${projectDir}/assets/consolidate_design_metrics.py", checkIfExists: true).first()

        // Collect output files from each analysis process
        // These will be staged into the consolidation task's work directory

        // ipSAE scores (the .txt files, not byres)
        ch_ipsae_files = (params.run_ipsae && params.run_proteinmpnn && params.run_boltz2_refold)
            ? IPSAE_CALCULATE.out.scores
                .map { meta, file -> file }
                .collect()
                .ifEmpty { file('NO_IPSAE_FILES') }
            : Channel.value(file('NO_IPSAE_FILES'))

        // Prodigy results (.txt files)
        ch_prodigy_files = (params.run_prodigy && params.run_proteinmpnn && params.run_boltz2_refold)
            ? PRODIGY_PREDICT.out.results
                .map { meta, file -> file }
                .collect()
                .ifEmpty { file('NO_PRODIGY_FILES') }
            : Channel.value(file('NO_PRODIGY_FILES'))

        // Foldseek summaries (.tsv files)
        ch_foldseek_files = (params.run_foldseek && params.run_proteinmpnn && params.run_boltz2_refold)
            ? FOLDSEEK_SEARCH.out.summary
                .map { meta, file -> file }
                .collect()
                .ifEmpty { file('NO_FOLDSEEK_FILES') }
            : Channel.value(file('NO_FOLDSEEK_FILES'))

        // Run consolidation with staged files
        CONSOLIDATE_METRICS(
            ch_ipsae_files,
            ch_prodigy_files,
            ch_foldseek_files,
            ch_consolidate_script
        )
    }

    emit:
    // Boltzgen outputs
    boltzgen_results = BOLTZGEN_RUN.out.results
    final_designs = BOLTZGEN_RUN.out.final_designs
    
    // ProteinMPNN outputs (will be empty if not run)
    mpnn_optimized = params.run_proteinmpnn ? PROTEINMPNN_OPTIMIZE.out.optimized_designs : Channel.empty()
    mpnn_sequences = params.run_proteinmpnn ? PROTEINMPNN_OPTIMIZE.out.sequences : Channel.empty()
    mpnn_scores = params.run_proteinmpnn ? PROTEINMPNN_OPTIMIZE.out.scores : Channel.empty()
    
    // Boltz-2 refolding outputs (will be empty if not run)
    boltz2_structures = (params.run_proteinmpnn && params.run_boltz2_refold) ? BOLTZ2_REFOLD.out.structures : Channel.empty()
    boltz2_confidence = (params.run_proteinmpnn && params.run_boltz2_refold) ? BOLTZ2_REFOLD.out.confidence : Channel.empty()
    boltz2_pae_npz = (params.run_proteinmpnn && params.run_boltz2_refold) ? BOLTZ2_REFOLD.out.pae_npz : Channel.empty()
    boltz2_affinity = (params.run_proteinmpnn && params.run_boltz2_refold) ? BOLTZ2_REFOLD.out.affinity : Channel.empty()
    
    // Optional analysis outputs (will be empty if not run)
    foldseek_results = (params.run_foldseek && params.run_proteinmpnn && params.run_boltz2_refold) ? FOLDSEEK_SEARCH.out.results : Channel.empty()
    foldseek_summary = (params.run_foldseek && params.run_proteinmpnn && params.run_boltz2_refold) ? FOLDSEEK_SEARCH.out.summary : Channel.empty()
    
    // Consolidation outputs (will be empty if not run)
    metrics_summary = params.run_consolidation ? CONSOLIDATE_METRICS.out.summary_csv : Channel.empty()
    metrics_report = params.run_consolidation ? CONSOLIDATE_METRICS.out.report_html : Channel.empty()
}
