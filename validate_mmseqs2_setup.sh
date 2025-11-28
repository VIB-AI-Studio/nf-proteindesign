#!/bin/bash
################################################################################
# MMSeqs2 MSA Setup Validation Script
################################################################################
# This script validates the MMSeqs2 and GPU setup for the nf-proteindesign pipeline
# Run this before using MSA features to ensure everything is configured correctly
################################################################################

set -e

echo "============================================"
echo "MMSeqs2 MSA Setup Validation"
echo "============================================"
echo ""

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

ERRORS=0
WARNINGS=0

# Function to print status
print_status() {
    local status=$1
    local message=$2
    
    if [ "$status" = "OK" ]; then
        echo -e "${GREEN}✓${NC} $message"
    elif [ "$status" = "WARN" ]; then
        echo -e "${YELLOW}⚠${NC} $message"
        ((WARNINGS++))
    else
        echo -e "${RED}✗${NC} $message"
        ((ERRORS++))
    fi
}

################################################################################
# 1. Check MMSeqs2 Installation
################################################################################
echo "1. Checking MMSeqs2 installation..."
echo "-----------------------------------"

if command -v mmseqs &> /dev/null; then
    MMSEQS_VERSION=$(mmseqs version 2>&1 | head -n1 | awk '{print $2}')
    print_status "OK" "MMSeqs2 is installed: version $MMSEQS_VERSION"
    
    # Check version compatibility (need >= 13.45111 for GPU support)
    REQUIRED_VERSION="13.45111"
    if [[ "$(printf '%s\n' "$REQUIRED_VERSION" "$MMSEQS_VERSION" | sort -V | head -n1)" == "$REQUIRED_VERSION" ]]; then
        print_status "OK" "MMSeqs2 version is compatible (>= $REQUIRED_VERSION)"
    else
        print_status "WARN" "MMSeqs2 version might not support GPU (found: $MMSEQS_VERSION, recommended: >= $REQUIRED_VERSION)"
    fi
else
    print_status "ERROR" "MMSeqs2 is not installed or not in PATH"
    echo "         Install with: conda install -c bioconda mmseqs2"
fi

echo ""

################################################################################
# 2. Check GPU Availability
################################################################################
echo "2. Checking GPU availability..."
echo "-------------------------------"

if command -v nvidia-smi &> /dev/null; then
    if nvidia-smi &> /dev/null; then
        print_status "OK" "NVIDIA GPU detected"
        
        # Get GPU information
        GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -n1)
        GPU_MEMORY=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader | head -n1)
        CUDA_VERSION=$(nvidia-smi | grep "CUDA Version" | awk '{print $9}')
        
        echo "         GPU: $GPU_NAME"
        echo "         Memory: $GPU_MEMORY"
        echo "         CUDA: $CUDA_VERSION"
        
        # Check if GPU has enough memory (recommend >= 8GB)
        GPU_MEMORY_GB=$(echo $GPU_MEMORY | awk '{print $1}' | cut -d'.' -f1)
        if [ "$GPU_MEMORY_GB" -ge 8 ]; then
            print_status "OK" "GPU memory is sufficient (>= 8GB)"
        else
            print_status "WARN" "GPU memory is limited (< 8GB) - may fail on large sequences"
        fi
    else
        print_status "WARN" "NVIDIA driver detected but GPU not accessible"
        echo "         Check driver installation and permissions"
    fi
else
    print_status "WARN" "No NVIDIA GPU detected - will use CPU mode (slower)"
    echo "         MMSeqs2 MSA will still work but take ~10x longer"
fi

echo ""

################################################################################
# 3. Check Database Configuration
################################################################################
echo "3. Checking database configuration..."
echo "-------------------------------------"

# Check if database is specified in environment or config
if [ -n "$MMSEQS2_DB" ]; then
    DB_PATH="$MMSEQS2_DB"
    echo "Using database from MMSEQS2_DB environment variable"
elif [ -f "nextflow.config" ]; then
    # Try to extract from config file
    DB_PATH=$(grep "mmseqs2_database" nextflow.config | grep -v "//" | awk -F"=" '{print $2}' | tr -d ' "' | tr -d "'")
    if [ -n "$DB_PATH" ] && [ "$DB_PATH" != "null" ]; then
        echo "Using database from nextflow.config: $DB_PATH"
    else
        DB_PATH=""
    fi
else
    DB_PATH=""
fi

if [ -z "$DB_PATH" ]; then
    print_status "WARN" "MMSeqs2 database not configured"
    echo "         Set params.mmseqs2_database in nextflow.config or export MMSEQS2_DB"
    echo ""
    echo "         Download databases from:"
    echo "         - UniRef30 (~90GB): https://wwwuser.gwdg.de/~compbiol/colabfold/uniref30_2202_db.tar.gz"
    echo "         - ColabFoldDB (~1.5TB): https://wwwuser.gwdg.de/~compbiol/colabfold/colabfold_envdb_202108.tar.gz"
else
    if [ -d "$DB_PATH" ]; then
        print_status "OK" "Database directory found: $DB_PATH"
        
        # Check for database files
        if [ -f "${DB_PATH}.dbtype" ] || [ -f "${DB_PATH}/db.dbtype" ]; then
            print_status "OK" "Database appears to be properly formatted"
            
            # Get database size
            DB_SIZE=$(du -sh "$DB_PATH" 2>/dev/null | awk '{print $1}')
            echo "         Database size: $DB_SIZE"
            
            # Check if database is indexed
            if [ -f "${DB_PATH}.index" ] || [ -f "${DB_PATH}/db.index" ]; then
                print_status "OK" "Database is indexed (ready for fast searches)"
            else
                print_status "WARN" "Database not indexed - first search will be slower"
                echo "         Run: mmseqs createindex $DB_PATH tmp"
            fi
        else
            print_status "ERROR" "Database directory exists but appears incomplete"
            echo "         Missing .dbtype file - database may be corrupted"
        fi
    else
        print_status "ERROR" "Database path not found: $DB_PATH"
        echo "         Download and extract a database to this location"
    fi
fi

echo ""

################################################################################
# 4. Check Docker/Singularity GPU Support
################################################################################
echo "4. Checking container GPU support..."
echo "-------------------------------------"

# Check Docker
if command -v docker &> /dev/null; then
    print_status "OK" "Docker is installed"
    
    # Test GPU access
    if docker run --rm --gpus all nvidia/cuda:11.8.0-base-ubuntu22.04 nvidia-smi &> /dev/null; then
        print_status "OK" "Docker has GPU access"
    else
        print_status "WARN" "Docker cannot access GPU"
        echo "         Install nvidia-docker2: https://github.com/NVIDIA/nvidia-docker"
    fi
else
    print_status "WARN" "Docker not installed"
fi

# Check Singularity
if command -v singularity &> /dev/null; then
    print_status "OK" "Singularity is installed"
    SINGULARITY_VERSION=$(singularity --version)
    echo "         Version: $SINGULARITY_VERSION"
else
    print_status "WARN" "Singularity not installed (optional)"
fi

echo ""

################################################################################
# 5. Check Nextflow Configuration
################################################################################
echo "5. Checking Nextflow configuration..."
echo "--------------------------------------"

if command -v nextflow &> /dev/null; then
    NF_VERSION=$(nextflow -version 2>&1 | grep version | awk '{print $3}')
    print_status "OK" "Nextflow is installed: version $NF_VERSION"
    
    # Check if version is recent enough (need >= 22.10 for better GPU support)
    REQUIRED_NF="22.10"
    if [[ "$(printf '%s\n' "$REQUIRED_NF" "$NF_VERSION" | sort -V | head -n1)" == "$REQUIRED_NF" ]]; then
        print_status "OK" "Nextflow version is compatible (>= $REQUIRED_NF)"
    else
        print_status "WARN" "Nextflow version is old (found: $NF_VERSION, recommended: >= $REQUIRED_NF)"
        echo "         Update with: nextflow self-update"
    fi
else
    print_status "ERROR" "Nextflow is not installed"
    echo "         Install from: https://www.nextflow.io/docs/latest/getstarted.html"
fi

# Check config file
if [ -f "nextflow.config" ]; then
    print_status "OK" "nextflow.config found"
    
    # Check if MSA parameters are configured
    if grep -q "boltz2_msa_mode" nextflow.config; then
        MSA_MODE=$(grep "boltz2_msa_mode" nextflow.config | grep -v "//" | awk -F"=" '{print $2}' | tr -d ' "' | tr -d "'")
        print_status "OK" "MSA mode configured: $MSA_MODE"
    else
        print_status "WARN" "MSA parameters not found in config (will use defaults)"
    fi
else
    print_status "WARN" "nextflow.config not found in current directory"
fi

echo ""

################################################################################
# 6. System Resources Check
################################################################################
echo "6. Checking system resources..."
echo "-------------------------------"

# Check available memory
TOTAL_MEM=$(free -g | awk '/^Mem:/ {print $2}')
AVAIL_MEM=$(free -g | awk '/^Mem:/ {print $7}')
print_status "OK" "Total memory: ${TOTAL_MEM}GB, Available: ${AVAIL_MEM}GB"

if [ "$AVAIL_MEM" -lt 16 ]; then
    print_status "WARN" "Available memory is low (< 16GB) - may need to reduce batch sizes"
fi

# Check CPU cores
CPU_CORES=$(nproc)
print_status "OK" "CPU cores: $CPU_CORES"

# Check disk space
DISK_FREE=$(df -h . | awk 'NR==2 {print $4}')
print_status "OK" "Available disk space: $DISK_FREE"

echo ""

################################################################################
# Summary
################################################################################
echo "============================================"
echo "Validation Summary"
echo "============================================"

if [ $ERRORS -eq 0 ] && [ $WARNINGS -eq 0 ]; then
    echo -e "${GREEN}✓ All checks passed!${NC} System is ready for MMSeqs2 MSA generation."
elif [ $ERRORS -eq 0 ]; then
    echo -e "${YELLOW}⚠ Validation completed with $WARNINGS warning(s).${NC}"
    echo "  The system will work but with reduced performance or features."
else
    echo -e "${RED}✗ Validation failed with $ERRORS error(s) and $WARNINGS warning(s).${NC}"
    echo "  Fix the errors above before running the pipeline."
fi

echo ""
echo "Quick Start Command:"
echo "-------------------"
echo "nextflow run main.nf \\"
echo "  --input samplesheet.csv \\"
echo "  --outdir results \\"
echo "  --run_proteinmpnn true \\"
echo "  --run_boltz2_refold true \\"
echo "  --boltz2_msa_mode target_only \\"
echo "  --mmseqs2_database /path/to/uniref30_2202_db"
echo ""
echo "For more information, see: docs/mmseqs2_msa_implementation.md"
echo "============================================"

exit $ERRORS
