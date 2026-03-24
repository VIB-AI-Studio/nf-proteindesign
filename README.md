# 🧬 nf-proteindesign
> ⚠️ **IMPORTANT**: This pipeline was developed by Seqera as a proof of principle using Seqera AI. It demonstrates the capabilities of AI-assisted bioinformatics pipeline development but should be thoroughly validated before use in production environments.

A Nextflow pipeline for AI-powered protein design using Boltzgen to design protein binders, nanobodies, and peptides.

## 📋 Overview
This pipeline automates the process of designing novel protein binders using Boltzgen and provides comprehensive analysis through optional modules:
- 🎯 **Boltzgen Design**: Generate protein, nanobody, or peptide binders for target structures
- 🧬 **ProteinMPNN**: Optimize sequences for improved stability and expression
- 🔄 **Boltz-2 Refolding**: Validate designs through structure prediction
- 📊 **IPSAE**: Score protein-protein interface quality
- ⚡ **PRODIGY**: Predict binding affinity
- 🔍 **Foldseek**: Search structural databases for similar designs
- 📈 **Metrics Consolidation**: Generate comprehensive analysis reports

## 🚀 Quick Start

### ✅ Prerequisites
- ⚙️ Nextflow (≥23.10)
- 🐳 Docker, Singularity, or Apptainer
- 🎮 **GPU required** for optimal performance (H100, A100, L40S, etc.)
- 💾 Sufficient storage for model weights (~6GB for Boltzgen, ~6GB for Boltz2)

### 🧪 Running with Test Profiles

#### Local/Single Node
Test the pipeline with one of three available profiles:
```bash
# Test protein binder design (recommended starting point)
nextflow run main.nf -profile test_design_protein,apptainer

# Test nanobody binder design  
nextflow run main.nf -profile test_design_nanobody,apptainer

# Test peptide binder design
nextflow run main.nf -profile test_design_peptide,apptainer
```

#### SLURM Cluster (Recommended for Production)
For multi-node clusters with GPU partitions:
```bash
# Flexible multi-queue scheduling (SLURM automatically picks best nodes)
nextflow run main.nf -profile test_design_protein,apptainer,slurm_flexible \
  --cache_dir /path/to/shared/cache

# Monitor jobs with: squeue -u $(whoami)
```

### 🔧 Container & Configuration Notes

#### Container Timeouts
Large containers may require extended pull timeouts:
```bash
# For Singularity/Apptainer (increase as needed)
nextflow run main.nf -profile test_design_protein,apptainer \
  --apptainer.pullTimeout '120m' \
  --cache_dir /path/to/cache
```

#### Boltz-2 Affinity Prediction Limitations
**Important**: Boltz-2's affinity prediction only works with small molecule ligands, NOT protein-protein interactions. The test profiles have been configured to disable affinity prediction by default since protein binders, nanobodies, and peptides are all treated as proteins, not ligands by Boltz-2.

- ✅ **Structure prediction**: Works for all design types
- ❌ **Affinity prediction**: Only for small molecules (disabled in test profiles)  
- ✅ **Alternative**: Use PRODIGY for protein-protein binding affinity predictions

### 🔬 Running with Your Own Data
```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results \
  -profile apptainer,slurm_flexible
```

## 📝 Input Format
The pipeline requires a CSV samplesheet with design specifications. See `assets/test_data/` for examples:
```csv
sample,design_yaml,protocol,num_designs,budget
my_design,design.yaml,protein-anything,10,5
```

## ⚙️ Key Parameters

### Core Pipeline Parameters
- `--input`: Path to samplesheet CSV
- `--outdir`: Output directory (default: `./results`)
- `--cache_dir`: Cache directory for model weights (recommended: shared filesystem location, ~12GB total)

### Analysis Module Controls
- `--run_proteinmpnn`: Enable ProteinMPNN sequence optimization (default: `true`)
- `--run_boltz2_refold`: Enable Boltz-2 structure prediction (default: `true`)  
- `--run_ipsae`: Enable IPSAE interface scoring (default: `true`)
- `--run_prodigy`: Enable PRODIGY affinity prediction (default: `true`)
- `--run_consolidation`: Generate consolidated metrics report (default: `true`)

### Boltz-2 Specific Parameters
- `--boltz2_predict_affinity`: Enable affinity prediction (default: `true`, but disabled in test profiles due to ligand limitations)
- `--boltz2_use_msa`: Use MSA server for designs without MSA files (default: `false`, enabled for peptides)
- `--boltz2_num_recycling`: Number of recycling steps (default: `3`)
- `--boltz2_num_diffusion`: Number of diffusion samples (default: `200`)

### Container & HPC Parameters  
- `--apptainer.pullTimeout`: Container pull timeout (default: `20m`, recommend `120m+` for slow networks)
- `--singularity.pullTimeout`: Alternative for Singularity (default: `20m`)

See `nextflow.config` for all available parameters.

## 🏗️ Infrastructure Support

### 🖥️ Local/Single Node
```bash
# Docker (if available)
nextflow run main.nf -profile test_design_protein,docker

# Apptainer/Singularity (recommended for HPC)
nextflow run main.nf -profile test_design_protein,apptainer
```

### 🚀 SLURM Clusters  
The pipeline includes a flexible SLURM profile (`slurm_flexible`) that:
- ✅ **Auto-distributes** GPU jobs across available GPU partitions
- ✅ **Load-balances** CPU jobs across available CPU nodes  
- ✅ **Optimizes costs** by using appropriate node types for each task
- ✅ **Handles failures** gracefully with multi-queue submission

```bash
# SLURM with flexible queue scheduling
nextflow run main.nf -profile test_design_protein,apptainer,slurm_flexible
```

### 🎮 GPU Requirements
- **Required for**: Boltzgen design, Boltz-2 refolding
- **Optional for**: ProteinMPNN optimization, Foldseek searches  
- **Recommended**: H100, A100, L40S, or similar NVIDIA GPUs
- **Memory**: 16GB+ GPU memory recommended for large designs

## 🔄 Resume & Recovery  
Nextflow automatically supports resume functionality:
```bash
# Resume from last checkpoint (same parameters required)
nextflow run main.nf -profile test_design_protein,apptainer,slurm_flexible \
  -resume

# Resume from specific run hash
nextflow run main.nf -profile test_design_protein,apptainer,slurm_flexible \
  -resume [RUN_HASH]

# Check previous runs
nextflow log
```

## 🐛 Troubleshooting

### Common Issues & Solutions

1. **Container timeout**: Increase `--apptainer.pullTimeout` to `120m` or higher for slow networks
2. **GPU not detected**: Ensure `--nv` flag and GPU drivers are properly configured  
3. **Out of space during container pull**: 
   - Set `--cache_dir` to location with sufficient storage (~12GB for model weights)
   - Or redirect Apptainer cache: `export APPTAINER_CACHEDIR=/path/to/large/storage`
   - Default cache location: `~/.apptainer/cache/` (often limited space)
4. **SLURM job failures**: Check queue availability with `sinfo` and `squeue`
5. **Boltz-2 affinity errors**: This is expected - affinity prediction only works for small molecule ligands, not protein-protein interactions. The test profiles disable this by default.
6. **Missing MSA errors**: For peptides/designs without MSA files, ensure `boltz2_use_msa = true` in your configuration

### Debug Commands
```bash
# Check available SLURM partitions
sinfo

# Monitor running jobs
squeue -u $(whoami)

# Test GPU access in container  
apptainer exec --nv [CONTAINER] nvidia-smi

# Check Nextflow run history
nextflow log -f name,status,duration,hash

# Resume from specific failed run
nextflow run main.nf -profile [PROFILES] -resume [HASH]
```

### Test Profile Configuration Issues
If you encounter issues with the test profiles, the following parameters are pre-configured:
- `boltz2_predict_affinity = false` - Prevents ligand classification errors
- `boltz2_use_msa = true` - Enables MSA server for designs without MSA files (peptides)
- Appropriate resource allocation for different compute environments

## 📁 Output
Results are organized by sample in the output directory:
```
results/
├── boltzgen/          # Boltzgen designs and structures
├── proteinmpnn/       # Optimized sequences (if enabled)
├── boltz2/            # Refolded structures (if enabled)
├── ipsae/             # Interface scores (if enabled)
├── prodigy/           # Affinity predictions (if enabled)
├── foldseek/          # Structural search results (if enabled)
└── consolidated/      # Combined metrics report (if enabled)
```

## 📚 Citation
If you use this pipeline, please cite:
- **Boltzgen**: [Add Boltzgen citation]
- **ProteinMPNN**: Dauparas et al. (2022) Science  
- **Boltz-2**: [Add Boltz-2 citation]
- **Nextflow**: Di Tommaso et al. (2017) Nature Biotechnology

## 🤝 Contributing
Found a bug or have suggestions? Please open an [issue](https://github.com/seqeralabs/nf-proteindesign/issues) or submit a pull request.

## 📄 License
This pipeline is distributed under the MIT License. See LICENSE for details.

---
<div align="center">

**Built with ❤️ using Nextflow and Seqera AI**

[Documentation](https://seqeralabs.github.io/nf-proteindesign/) | [Issues](https://github.com/seqeralabs/nf-proteindesign/issues) | [Seqera](https://seqera.io)

</div>