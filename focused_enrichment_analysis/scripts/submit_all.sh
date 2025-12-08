#!/bin/bash
#SBATCH --job-name=submit_all
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/submit_all.out
#SBATCH --error=logs/submit_all.err


################################################################################
# Master submission script for comprehensive enrichment analysis
################################################################################

# Usage message
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Submit enrichment analysis jobs to SLURM

OPTIONS:
    -h, --help              Show this help message
    -a, --all               Submit all analyses (data loading + all approaches + tiers)
    -d, --data-only         Submit only data loading
    -A, --approaches-only   Submit only approaches (1-6)
    -T, --tiers-only        Submit only tiers (1-2)
    -1 to -6                Submit specific approach (e.g., -1 for approach1)
    -t1, -t2                Submit specific tier
    --sequential            Run jobs sequentially (with dependencies)
    --parallel              Run jobs in parallel (default for approaches)

EXAMPLES:
    # Submit all analyses sequentially
    $0 --all --sequential

    # Submit only data loading
    $0 --data-only

    # Submit specific approach
    $0 -2

    # Submit all approaches in parallel
    $0 --approaches-only

    # Submit tier 1 after all approaches complete
    $0 -t1
EOF
    exit 1
}

# Parse command line arguments
ALL=0
DATA_ONLY=0
APPROACHES_ONLY=0
TIERS_ONLY=0
SEQUENTIAL=0
APPROACH1=0
APPROACH2=0
APPROACH3=0
APPROACH4=0
APPROACH5=0
APPROACH6=0
TIER1=0
TIER2=0

if [ $# -eq 0 ]; then
    usage
fi

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            usage
            ;;
        -a|--all)
            ALL=1
            shift
            ;;
        -d|--data-only)
            DATA_ONLY=1
            shift
            ;;
        -A|--approaches-only)
            APPROACHES_ONLY=1
            shift
            ;;
        -T|--tiers-only)
            TIERS_ONLY=1
            shift
            ;;
        --sequential)
            SEQUENTIAL=1
            shift
            ;;
        --parallel)
            SEQUENTIAL=0
            shift
            ;;
        -1)
            APPROACH1=1
            shift
            ;;
        -2)
            APPROACH2=1
            shift
            ;;
        -3)
            APPROACH3=1
            shift
            ;;
        -4)
            APPROACH4=1
            shift
            ;;
        -5)
            APPROACH5=1
            shift
            ;;
        -6)
            APPROACH6=1
            shift
            ;;
        -t1)
            TIER1=1
            shift
            ;;
        -t2)
            TIER2=1
            shift
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Set working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top

# Store job IDs
LOAD_DATA_JOB=""
APPROACH_JOBS=()
TIER_JOBS=()

################################################################################
# SUBMIT JOBS
################################################################################

echo "========================================================================"
echo "Submitting Comprehensive Enrichment Analysis Jobs"
echo "========================================================================"
echo ""

# 1. Data loading (always first if running anything)
if [ $ALL -eq 1 ] || [ $DATA_ONLY -eq 1 ] || [ $APPROACHES_ONLY -eq 1 ]; then
    echo "Submitting data loading job..."
    LOAD_DATA_JOB=$(sbatch --parsable SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_00_load_data.sh)
    echo "  Job ID: $LOAD_DATA_JOB"
    echo ""
fi

# 2. Approaches
if [ $ALL -eq 1 ] || [ $APPROACHES_ONLY -eq 1 ] || [ $APPROACH1 -eq 1 ]; then
    echo "Submitting Approach 1..."
    if [ ! -z "$LOAD_DATA_JOB" ]; then
        JOB_ID=$(sbatch --parsable --dependency=afterok:$LOAD_DATA_JOB \
                 SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_approach1.sh)
    else
        JOB_ID=$(sbatch --parsable SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_approach1.sh)
    fi
    APPROACH_JOBS+=($JOB_ID)
    echo "  Job ID: $JOB_ID"
    echo ""
fi

if [ $ALL -eq 1 ] || [ $APPROACHES_ONLY -eq 1 ] || [ $APPROACH2 -eq 1 ]; then
    echo "Submitting Approach 2..."
    if [ $SEQUENTIAL -eq 1 ] && [ ${#APPROACH_JOBS[@]} -gt 0 ]; then
        JOB_ID=$(sbatch --parsable --dependency=afterok:${APPROACH_JOBS[-1]} \
                 SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_approach2.sh)
    elif [ ! -z "$LOAD_DATA_JOB" ]; then
        JOB_ID=$(sbatch --parsable --dependency=afterok:$LOAD_DATA_JOB \
                 SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_approach2.sh)
    else
        JOB_ID=$(sbatch --parsable SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_approach2.sh)
    fi
    APPROACH_JOBS+=($JOB_ID)
    echo "  Job ID: $JOB_ID"
    echo ""
fi

if [ $ALL -eq 1 ] || [ $APPROACHES_ONLY -eq 1 ] || [ $APPROACH3 -eq 1 ]; then
    echo "Submitting Approach 3..."
    if [ $SEQUENTIAL -eq 1 ] && [ ${#APPROACH_JOBS[@]} -gt 0 ]; then
        JOB_ID=$(sbatch --parsable --dependency=afterok:${APPROACH_JOBS[-1]} \
                 SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_approach3.sh)
    elif [ ! -z "$LOAD_DATA_JOB" ]; then
        JOB_ID=$(sbatch --parsable --dependency=afterok:$LOAD_DATA_JOB \
                 SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_approach3.sh)
    else
        JOB_ID=$(sbatch --parsable SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_approach3.sh)
    fi
    APPROACH_JOBS+=($JOB_ID)
    echo "  Job ID: $JOB_ID"
    echo ""
fi

if [ $ALL -eq 1 ] || [ $APPROACHES_ONLY -eq 1 ] || [ $APPROACH4 -eq 1 ]; then
    echo "Submitting Approach 4..."
    if [ $SEQUENTIAL -eq 1 ] && [ ${#APPROACH_JOBS[@]} -gt 0 ]; then
        JOB_ID=$(sbatch --parsable --dependency=afterok:${APPROACH_JOBS[-1]} \
                 SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_approach4.sh)
    elif [ ! -z "$LOAD_DATA_JOB" ]; then
        JOB_ID=$(sbatch --parsable --dependency=afterok:$LOAD_DATA_JOB \
                 SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_approach4.sh)
    else
        JOB_ID=$(sbatch --parsable SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_approach4.sh)
    fi
    APPROACH_JOBS+=($JOB_ID)
    echo "  Job ID: $JOB_ID"
    echo ""
fi

if [ $ALL -eq 1 ] || [ $APPROACHES_ONLY -eq 1 ] || [ $APPROACH5 -eq 1 ]; then
    echo "Submitting Approach 5..."
    if [ $SEQUENTIAL -eq 1 ] && [ ${#APPROACH_JOBS[@]} -gt 0 ]; then
        JOB_ID=$(sbatch --parsable --dependency=afterok:${APPROACH_JOBS[-1]} \
                 SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_approach5.sh)
    elif [ ! -z "$LOAD_DATA_JOB" ]; then
        JOB_ID=$(sbatch --parsable --dependency=afterok:$LOAD_DATA_JOB \
                 SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_approach5.sh)
    else
        JOB_ID=$(sbatch --parsable SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_approach5.sh)
    fi
    APPROACH_JOBS+=($JOB_ID)
    echo "  Job ID: $JOB_ID"
    echo ""
fi

if [ $ALL -eq 1 ] || [ $APPROACHES_ONLY -eq 1 ] || [ $APPROACH6 -eq 1 ]; then
    echo "Submitting Approach 6..."
    if [ $SEQUENTIAL -eq 1 ] && [ ${#APPROACH_JOBS[@]} -gt 0 ]; then
        JOB_ID=$(sbatch --parsable --dependency=afterok:${APPROACH_JOBS[-1]} \
                 SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_approach6.sh)
    elif [ ! -z "$LOAD_DATA_JOB" ]; then
        JOB_ID=$(sbatch --parsable --dependency=afterok:$LOAD_DATA_JOB \
                 SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_approach6.sh)
    else
        JOB_ID=$(sbatch --parsable SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_approach6.sh)
    fi
    APPROACH_JOBS+=($JOB_ID)
    echo "  Job ID: $JOB_ID"
    echo ""
fi

# 3. Tiers (run after all approaches if --all)
if [ $ALL -eq 1 ] || [ $TIERS_ONLY -eq 1 ] || [ $TIER1 -eq 1 ]; then
    echo "Submitting Tier 1..."
    if [ ${#APPROACH_JOBS[@]} -gt 0 ]; then
        # Create dependency string for all approach jobs
        DEPS=$(IFS=:; echo "${APPROACH_JOBS[*]}")
        JOB_ID=$(sbatch --parsable --dependency=afterok:$DEPS \
                 SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_tier1.sh)
    else
        JOB_ID=$(sbatch --parsable SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_tier1.sh)
    fi
    TIER_JOBS+=($JOB_ID)
    echo "  Job ID: $JOB_ID"
    echo ""
fi

if [ $ALL -eq 1 ] || [ $TIERS_ONLY -eq 1 ] || [ $TIER2 -eq 1 ]; then
    echo "Submitting Tier 2..."
    if [ ${#APPROACH_JOBS[@]} -gt 0 ]; then
        # Create dependency string for all approach jobs
        DEPS=$(IFS=:; echo "${APPROACH_JOBS[*]}")
        JOB_ID=$(sbatch --parsable --dependency=afterok:$DEPS \
                 SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_tier2.sh)
    else
        JOB_ID=$(sbatch --parsable SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/submit_tier2.sh)
    fi
    TIER_JOBS+=($JOB_ID)
    echo "  Job ID: $JOB_ID"
    echo ""
fi

################################################################################
# SUMMARY
################################################################################

echo "========================================================================"
echo "SUBMISSION COMPLETE"
echo "========================================================================"
echo ""
echo "Submitted jobs:"
if [ ! -z "$LOAD_DATA_JOB" ]; then
    echo "  Data loading: $LOAD_DATA_JOB"
fi
if [ ${#APPROACH_JOBS[@]} -gt 0 ]; then
    echo "  Approaches: ${APPROACH_JOBS[@]}"
fi
if [ ${#TIER_JOBS[@]} -gt 0 ]; then
    echo "  Tiers: ${TIER_JOBS[@]}"
fi
echo ""
echo "Monitor jobs with: squeue -u $USER"
echo "View logs in: SRF_Eva_integrated_analysis/focused_enrichment_analysis/logs/"
echo ""
echo "========================================================================"
