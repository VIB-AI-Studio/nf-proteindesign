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
include { EXTRACT_TARGET_SEQUENCES } from '../modules/local/extract_target_sequences'
include { MMSEQS2_MSA } from '../modules/local/mmseqs2_msa'
include { BOLTZ2_REFOLD } from '../modules/local/boltz2_refold'
include { IPSAE_CALCULATE } from '../modules/local/ipsae_calculate'
include { PRODIGY_PREDICT } from '../modules/local/prodigy_predict'
include { FOLDSEEK_SEARCH } from '../modules/local/foldseek_search'
include { CONSOLIDATE_METRICS } from '../modules/local/consolidate_metrics'

workflow PROTEIN_DESIGN {
    
    take:
    ch_input    // channel: [meta, design_yaml, structure_files]
    ch_cache    // channel: path to cache directory or EMPTY_CACHE placeholder

    main:
    
    // ========================================================================
    // Run Boltzgen on design YAMLs
    // ========================================================================
    
    // Run Boltzgen for each design in parallel
    BOLTZGEN_RUN(ch_input, ch_cache)
    
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
                    def design_meta = [
                        id: "${meta.id}_${pdb_file.baseName}",
                        parent_id: meta.id,
                        design_name: pdb_file.baseName
                    ]
                    
                    [design_meta, pdb_file]
                }
            }
        
        // Run ProteinMPNN on each design individually (parallel execution per budget design)
        PROTEINMPNN_OPTIMIZE(ch_pdb_per_design)
        
        // Use ProteinMPNN optimized structures for downstream analyses
        ch_final_designs_for_analysis = PROTEINMPNN_OPTIMIZE.out.optimized_designs
        
        // ====================================================================
        // Step 3: Extract target sequences for Protenix refolding
        // ====================================================================
        // PURPOSE: Extract the TARGET sequence (binding partner) from Boltzgen structures
        // WHY: Protenix needs to know which chain is the target (to keep fixed) when
        //      refolding ProteinMPNN-optimized binder sequences
        // WHAT: Reads original Boltzgen CIF files and extracts the target chain sequence
        // OUTPUT: Plain text file with target sequence (one per design)
        // SOLUTION: Use ONLY the first budget design CIF (rank_1) to avoid naming collisions
        //           since the target sequence is identical across all designs for a sample
        // ====================================================================
        if (params.run_boltz2_refold) {
            // Extract target sequences from the FIRST budget design only
            // The target is the same across all designs, so we only need to extract it once
            ch_boltzgen_structures = BOLTZGEN_RUN.out.budget_design_cifs
                .map { meta, cif_files ->
                    // Take only the FIRST CIF file (typically rank_1.cif)
                    def cif_list = cif_files instanceof List ? cif_files : [cif_files]
                    def first_cif = cif_list.sort()[0]  // Sort to ensure consistent selection
                    [meta, first_cif]
                }
            
            EXTRACT_TARGET_SEQUENCES(ch_boltzgen_structures)
            
            // ================================================================
            // MMSeqs2 MSA Generation (if enabled)
            // ================================================================
            def msa_mode = params.boltz2_msa_mode ?: 'target_only'
            
            // Prepare target MSAs if needed
            if (msa_mode in ['target_only', 'both']) {
                // Create target sequence FASTA files for MMSeqs2
                // Deduplicate targets to run MSA only once per unique sequence
                ch_target_for_msa = EXTRACT_TARGET_SEQUENCES.out.target_sequences
                    .map { meta, target_seq ->
                        // Read sequence content for deduplication
                        def seq_content = target_seq.text.trim()
                        [seq_content, meta, target_seq]
                    }
                    .unique { it[0] }  // Deduplicate by sequence content
                    .map { seq_content, meta, target_seq ->
                        // Create meta for MSA process
                        def msa_meta = [
                            id: "${meta.id}_target",
                            parent_id: meta.id,
                            sequence_type: 'target'
                        ]
                        [msa_meta, target_seq]
                    }
                
                // Run MMSeqs2 on unique target sequences
                MMSEQS2_MSA(ch_target_for_msa)
                ch_target_msa = MMSEQS2_MSA.out.msa
            } else {
                // Create empty channel with NO_MSA placeholder
                ch_target_msa = Channel.empty()
            }
            
            // Prepare binder MSAs if needed
            if (msa_mode in ['binder_only', 'both']) {
                // Extract binder sequences from FASTA files and deduplicate
                ch_binder_for_msa = PROTEINMPNN_OPTIMIZE.out.sequences
                    .flatMap { meta, fasta_files ->
                        def fasta_list = fasta_files instanceof List ? fasta_files : [fasta_files]
                        fasta_list.collect { fasta_file ->
                            [meta, fasta_file]
                        }
                    }
                    .map { meta, fasta_file ->
                        // Read and parse FASTA to extract first sequence
                        def lines = fasta_file.readLines()
                        def seq_lines = []
                        def in_seq = false
                        for (line in lines) {
                            if (line.startsWith('>')) {
                                if (in_seq) break  // Stop after first sequence
                                in_seq = true
                            } else if (in_seq) {
                                seq_lines.add(line.trim())
                            }
                        }
                        def binder_seq = seq_lines.join('')
                        [binder_seq, meta, fasta_file]
                    }
                    .unique { it[0] }  // Deduplicate by sequence content
                    .map { binder_seq, meta, fasta_file ->
                        // Create temporary FASTA for MSA
                        def msa_meta = [
                            id: "${meta.id}_binder",
                            parent_id: meta.parent_id,
                            sequence_type: 'binder'
                        ]
                        
                        // Create a proper FASTA file from the sequence
                        def temp_fasta = file("${workDir}/binder_seqs/${meta.id}_binder.fasta")
                        temp_fasta.getParent().mkdirs()
                        temp_fasta.text = ">${meta.id}_binder\\n${binder_seq}\\n"
                        
                        [msa_meta, temp_fasta]
                    }
                
                // Run MMSeqs2 on unique binder sequences
                MMSEQS2_MSA(ch_binder_for_msa)
                ch_binder_msa = MMSEQS2_MSA.out.msa
            } else {
                // Create empty channel with NO_MSA placeholder
                ch_binder_msa = Channel.empty()
            }
            
            // ================================================================
            // Prepare inputs for Boltz-2 with MSA support
            // ================================================================
            // Parallelize Boltz-2 per FASTA file (one per ProteinMPNN sequence)
            ch_boltz2_base = PROTEINMPNN_OPTIMIZE.out.sequences
                .flatMap { meta, fasta_files ->
                    def fasta_list = fasta_files instanceof List ? new ArrayList(fasta_files) : [fasta_files]
                    fasta_list.collect { fasta_file ->
                        def seq_meta = [
                            id: "${meta.id}_${fasta_file.baseName}",
                            parent_id: meta.parent_id,
                            mpnn_parent_id: meta.id,
                            sequence_name: fasta_file.baseName
                        ]
                        [seq_meta, fasta_file]
                    }
                }
                .map { meta, fasta -> 
                    [meta.parent_id, meta, fasta]
                }
                .combine(
                    EXTRACT_TARGET_SEQUENCES.out.target_sequences.map { meta, seq -> 
                        [meta.id, seq]
                    },
                    by: 0
                )
                .map { parent_id, meta, fasta, target_seq ->
                    [meta, fasta, target_seq]
                }
            
            // Add MSA files based on mode
            if (msa_mode == 'target_only') {
                ch_boltz2_input = ch_boltz2_base
                    .map { meta, fasta, target_seq ->
                        [meta.parent_id, meta, fasta, target_seq]
                    }
                    .combine(
                        ch_target_msa.map { msa_meta, msa_file ->
                            // Extract parent ID (remove _target suffix)
                            def parent_id = msa_meta.id.replaceAll(/_target$/, '')
                            [parent_id, msa_file]
                        },
                        by: 0
                    )
                    .map { parent_id, meta, fasta, target_seq, target_msa ->
                        [meta, fasta, target_seq, target_msa, file('NO_MSA')]
                    }
            } else if (msa_mode == 'binder_only') {
                ch_boltz2_input = ch_boltz2_base
                    .map { meta, fasta, target_seq ->
                        [meta.mpnn_parent_id, meta, fasta, target_seq]
                    }
                    .combine(
                        ch_binder_msa.map { msa_meta, msa_file ->
                            def parent_id = msa_meta.id.replaceAll(/_binder$/, '')
                            [parent_id, msa_file]
                        },
                        by: 0
                    )
                    .map { mpnn_parent_id, meta, fasta, target_seq, binder_msa ->
                        [meta, fasta, target_seq, file('NO_MSA'), binder_msa]
                    }
            } else if (msa_mode == 'both') {
                ch_boltz2_input = ch_boltz2_base
                    .map { meta, fasta, target_seq ->
                        [meta.parent_id, meta.mpnn_parent_id, meta, fasta, target_seq]
                    }
                    .combine(
                        ch_target_msa.map { msa_meta, msa_file ->
                            def parent_id = msa_meta.id.replaceAll(/_target$/, '')
                            [parent_id, msa_file]
                        },
                        by: 0
                    )
                    .map { parent_id, mpnn_parent_id, meta, fasta, target_seq, target_msa ->
                        [mpnn_parent_id, meta, fasta, target_seq, target_msa]
                    }
                    .combine(
                        ch_binder_msa.map { msa_meta, msa_file ->
                            def parent_id = msa_meta.id.replaceAll(/_binder$/, '')
                            [parent_id, msa_file]
                        },
                        by: 0
                    )
                    .map { mpnn_parent_id, meta, fasta, target_seq, target_msa, binder_msa ->
                        [meta, fasta, target_seq, target_msa, binder_msa]
                    }
            } else {
                // No MSA mode
                ch_boltz2_input = ch_boltz2_base
                    .map { meta, fasta, target_seq ->
                        [meta, fasta, target_seq, file('NO_MSA'), file('NO_MSA')]
                    }
            }
            
            // Run Boltz-2 structure prediction with MSA support
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
        // Prepare IPSAE script as a channel
        ch_ipsae_script = Channel.fromPath("${projectDir}/assets/ipsae.py", checkIfExists: true)
        
        // ====================================================================
        // Part 1: Process Boltzgen budget design CIF and NPZ files
        // ====================================================================
        // Process ALL budget design CIF and NPZ files from intermediate_designs_inverse_folded
        // This ensures we run IPSAE on ALL designs before filtering (e.g., if budget=10, run 10 times)
        ch_ipsae_boltzgen = BOLTZGEN_RUN.out.budget_design_cifs
            .join(BOLTZGEN_RUN.out.budget_design_npz, by: 0)
            .flatMap { meta, cif_files, npz_files ->
                // Convert to list if single file and create defensive copies
                def cif_list = cif_files instanceof List ? new ArrayList(cif_files) : [cif_files]
                def npz_list = npz_files instanceof List ? new ArrayList(npz_files) : [npz_files]
                
                // Create a map of basenames to files for quick lookup
                def npz_map = [:]
                npz_list.each { npz_file ->
                    npz_map[npz_file.baseName] = npz_file
                }
                
                // Match CIF files with corresponding NPZ files
                cif_list.collect { cif_file ->
                    def base_name = cif_file.baseName
                    def npz_file = npz_map[base_name]
                    
                    if (npz_file) {
                        def model_meta = [
                            id: "${meta.id}_${base_name}",
                            parent_id: meta.id,
                            model_id: "${meta.id}_${base_name}",
                            source: "boltzgen"
                        ]
                        
                        [model_meta, npz_file, cif_file]
                    } else {
                        log.warn "⚠️  No matching NPZ file found for ${cif_file.name} in design ${meta.id}"
                        null
                    }
                }.findAll { it != null }  // Remove null entries where no match was found
            }
        
        // ====================================================================
        // Part 2: Process Boltz-2 refolded structures (if enabled)
        // ====================================================================
        // Add Boltz-2 NPZ files if Boltz-2 refolding is enabled
        if (params.run_proteinmpnn && params.run_boltz2_refold) {
            // Get NPZ and CIF pairs directly from Boltz-2 (native NPZ output!)
            ch_ipsae_boltz2 = BOLTZ2_REFOLD.out.structures
                .join(BOLTZ2_REFOLD.out.pae_npz, by: 0)
                .flatMap { meta, cif_files, npz_files ->
                    // Convert to lists if single files
                    def cif_list = cif_files instanceof List ? cif_files : [cif_files]
                    def npz_list = npz_files instanceof List ? npz_files : [npz_files]
                    
                    // Create a map of basenames for matching
                    def npz_map = [:]
                    npz_list.each { npz_file ->
                        // Extract base name (without _pae suffix if present)
                        def base_name = npz_file.baseName.replaceAll(/_pae.*$/, '')
                        npz_map[base_name] = npz_file
                    }
                    
                    // Match CIF files with their NPZ files
                    cif_list.collect { cif_file ->
                        def base_name = cif_file.baseName
                        def npz_file = npz_map[base_name]
                        
                        if (npz_file) {
                            def ipsae_meta = [
                                id: "${meta.id}_${base_name}",
                                parent_id: meta.parent_id,
                                mpnn_parent_id: meta.mpnn_parent_id,
                                model_id: base_name,
                                source: "boltz2"
                            ]
                            
                            [ipsae_meta, npz_file, cif_file]
                        } else {
                            log.warn "⚠️  No matching NPZ file found for ${cif_file.name}"
                            null
                        }
                    }.findAll { it != null }
                }
            
            // Combine both Boltzgen and Boltz-2 inputs
            ch_ipsae_input = ch_ipsae_boltzgen.mix(ch_ipsae_boltz2)
        } else {
            // Only Boltzgen inputs
            ch_ipsae_input = ch_ipsae_boltzgen
        }
        
        // Run IPSAE calculation for all CIF/NPZ pairs (Boltzgen + Protenix)
        IPSAE_CALCULATE(ch_ipsae_input, ch_ipsae_script)
    }
    
    // ========================================================================
    // OPTIONAL: PRODIGY binding affinity prediction if enabled
    // ========================================================================
    if (params.run_prodigy) {
        // Prepare PRODIGY parser script as a channel
        ch_prodigy_script = Channel.fromPath("${projectDir}/assets/parse_prodigy_output.py", checkIfExists: true)
        
        // Use ALL budget design CIF files from intermediate_designs_inverse_folded
        // This ensures we run PRODIGY on ALL designs before filtering (e.g., if budget=10, run 10 times)
        // Strategy: Use flatMap to create individual tasks for each CIF file
        ch_prodigy_input = BOLTZGEN_RUN.out.budget_design_cifs
            .flatMap { meta, cif_files ->
                // Convert to list if single file and create defensive copy
                def cif_list = cif_files instanceof List ? new ArrayList(cif_files) : [cif_files]
                
                // Create a separate entry for each CIF file
                cif_list.collect { cif_file ->
                    def base_name = cif_file.baseName
                    // Create new meta map explicitly to avoid concurrent modification
                    def design_meta = [
                        id: "${meta.id}_${base_name}",
                        parent_id: meta.id,
                        source: "boltzgen"
                    ]
                    
                    [design_meta, cif_file]
                }
            }
        
        // Add Boltz-2-refolded structures if available
        if (params.run_proteinmpnn && params.run_boltz2_refold) {
            ch_prodigy_boltz2 = BOLTZ2_REFOLD.out.structures
                .flatMap { meta, cif_files ->
                    // Convert to list if single file and create defensive copy
                    def cif_list = cif_files instanceof List ? new ArrayList(cif_files) : [cif_files]
                    
                    // Create a separate entry for each CIF file
                    cif_list.collect { cif_file ->
                        def base_name = cif_file.baseName
                        // Create new meta map explicitly to avoid concurrent modification
                        def design_meta = [
                            id: "${meta.id}_${base_name}_boltz2",
                            parent_id: meta.parent_id,  // Maintain link to original Boltzgen design
                            mpnn_parent_id: meta.id,     // Track which ProteinMPNN design this came from
                            source: "boltz2"
                        ]
                        
                        [design_meta, cif_file]
                    }
                }
            
            // Combine both sources
            ch_prodigy_input = ch_prodigy_input.mix(ch_prodigy_boltz2)
        }
        
        // Run PRODIGY binding affinity prediction for all CIF files (Boltzgen + Boltz-2)
        PRODIGY_PREDICT(ch_prodigy_input, ch_prodigy_script)
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
        // Part 1: Process Boltzgen budget design CIF files
        // ====================================================================
        // Use ALL budget design CIF files from intermediate_designs_inverse_folded
        // This searches for homologs of the original Boltzgen-designed structures
        ch_foldseek_boltzgen = BOLTZGEN_RUN.out.budget_design_cifs
            .flatMap { meta, cif_files ->
                // Convert to list if single file and create defensive copy
                def cif_list = cif_files instanceof List ? new ArrayList(cif_files) : [cif_files]
                
                // Create a separate entry for each CIF file
                cif_list.collect { cif_file ->
                    def base_name = cif_file.baseName
                    // Create new meta map explicitly to avoid concurrent modification
                    def design_meta = [
                        id: "${meta.id}_${base_name}",
                        parent_id: meta.id,
                        source: "boltzgen"
                    ]
                    
                    [design_meta, cif_file]
                }
            }
        
        // ====================================================================
        // Part 2: Add Boltz-2 refolded structures (if enabled)
        // ====================================================================
        // Search for homologs of the Boltz-2 refolded structures with MPNN sequences
        if (params.run_proteinmpnn && params.run_boltz2_refold) {
            ch_foldseek_boltz2 = BOLTZ2_REFOLD.out.structures
                .flatMap { meta, cif_files ->
                    // Convert to list if single file and create defensive copy
                    def cif_list = cif_files instanceof List ? new ArrayList(cif_files) : [cif_files]
                    
                    // Create a separate entry for each CIF file
                    cif_list.collect { cif_file ->
                        def base_name = cif_file.baseName
                        def design_meta = [
                            id: "${meta.id}_${base_name}_boltz2",
                            parent_id: meta.parent_id,  // Link to original Boltzgen design
                            mpnn_parent_id: meta.id,     // Track which ProteinMPNN design this came from
                            source: "boltz2"
                        ]
                        
                        [design_meta, cif_file]
                    }
                }
            
            // Combine both Boltzgen and Boltz-2 structures
            ch_foldseek_input = ch_foldseek_boltzgen.mix(ch_foldseek_boltz2)
        } else {
            // Only Boltzgen structures
            ch_foldseek_input = ch_foldseek_boltzgen
        }
        
        // Run Foldseek structural search for all CIF files (Boltzgen + Boltz-2)
        // This runs in parallel with IPSAE and PRODIGY analyses
        FOLDSEEK_SEARCH(ch_foldseek_input, ch_foldseek_database)
    }
    
    // ========================================================================
    // CONSOLIDATION: Generate comprehensive metrics report
    // ========================================================================
    if (params.run_consolidation) {
        // Prepare consolidation script as a channel
        ch_consolidate_script = Channel.fromPath("${projectDir}/assets/consolidate_design_metrics.py", checkIfExists: true)
        
        // Create a trigger channel that waits for all analyses to complete
        // Collect all output channels that need to complete before consolidation
        // Start with Boltzgen results (always runs)
        ch_trigger = BOLTZGEN_RUN.out.results.collect()
        
        // Mix in other outputs based on what's enabled
        if (params.run_proteinmpnn) {
            ch_trigger = ch_trigger.mix(PROTEINMPNN_OPTIMIZE.out.optimized_designs.collect())
        }
        
        if (params.run_ipsae) {
            ch_trigger = ch_trigger.mix(IPSAE_CALCULATE.out.scores.collect())
        }
        
        if (params.run_prodigy) {
            ch_trigger = ch_trigger.mix(PRODIGY_PREDICT.out.summary.collect())
        }
        
        if (params.run_foldseek) {
            ch_trigger = ch_trigger.mix(FOLDSEEK_SEARCH.out.summary.collect())
        }
        
        // After all outputs are collected, create a single trigger
        // and map it to the output directory path
        ch_outdir = ch_trigger
            .collect()
            .map { params.outdir }
        
        // Run consolidation after all analyses are complete
        CONSOLIDATE_METRICS(ch_outdir, ch_consolidate_script)
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
    foldseek_results = params.run_foldseek ? FOLDSEEK_SEARCH.out.results : Channel.empty()
    foldseek_summary = params.run_foldseek ? FOLDSEEK_SEARCH.out.summary : Channel.empty()
    
    // Consolidation outputs (will be empty if not run)
    metrics_summary = params.run_consolidation ? CONSOLIDATE_METRICS.out.summary_csv : Channel.empty()
    metrics_report = params.run_consolidation ? CONSOLIDATE_METRICS.out.report_markdown : Channel.empty()
}
