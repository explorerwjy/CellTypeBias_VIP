#!/bin/bash
# Rerun all DDD_61 pipelines after gene weight file fix
# Runs 3 jobs in parallel (each using 10 processes)
# Usage: bash scripts/rerun_DDD61.sh

set -euo pipefail

BASE="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP"
GW="${BASE}/dat/GeneWeights/DDD.top61.gw.bgmr.csv"
SPECMAT_CENTER="${BASE}/dat/ExpMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.mean_centered.csv"
SPECMAT_QN="${BASE}/dat/ExpMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.qn.csv"
LOGDIR="${BASE}/logs/rerun_DDD61"
mkdir -p "$LOGDIR"

# Track overall success
FAILED_JOBS=""

run_pipeline() {
    local TOPDIR="$1"
    local ATYPE="$2"
    local SPECMAT="$3"
    local STEP1_ARGS="$4"
    local JOB_ID="$5"

    local RESDIR="${BASE}/results/${TOPDIR}/${ATYPE}"
    local NULL_WEIGHTS="${RESDIR}/null_weights/DDD_61_random_geneweights.csv"
    local NULL_BIAS="${RESDIR}/null_bias/DDD_61_null_bias.csv"
    local BIAS_OUT="${RESDIR}/DDD_61_bias_addP.csv"
    local LOGFILE="${LOGDIR}/job_${JOB_ID}_${TOPDIR}_${ATYPE}.log"

    echo "[JOB ${JOB_ID}] Starting: ${TOPDIR}/${ATYPE}" | tee -a "$LOGFILE"
    echo "[JOB ${JOB_ID}] $(date)" | tee -a "$LOGFILE"

    # Ensure output directories exist
    mkdir -p "${RESDIR}/null_weights" "${RESDIR}/null_bias"

    # Step 1: Generate null gene weights
    echo "[JOB ${JOB_ID}] Step 1: Generating null gene weights..." | tee -a "$LOGFILE"
    python "${BASE}/scripts/script_generate_geneweights.py" \
        --WeightDF "$GW" \
        --SpecMat "$SPECMAT" \
        --n_sims 10000 \
        --n_processes 10 \
        ${STEP1_ARGS} \
        --outfile "$NULL_WEIGHTS" \
        >> "$LOGFILE" 2>&1

    if [ $? -ne 0 ]; then
        echo "[JOB ${JOB_ID}] FAILED at Step 1" | tee -a "$LOGFILE"
        return 1
    fi

    # Step 2: Compute null bias
    echo "[JOB ${JOB_ID}] Step 2: Computing null bias..." | tee -a "$LOGFILE"
    python "${BASE}/scripts/script_run_ctrl_sim.py" \
        --SpecMat "$SPECMAT" \
        --Ctrl_Genes_Fil "$NULL_WEIGHTS" \
        --outfile "$NULL_BIAS" \
        --n_processes 10 \
        -m human_ct_bias \
        >> "$LOGFILE" 2>&1

    if [ $? -ne 0 ]; then
        echo "[JOB ${JOB_ID}] FAILED at Step 2" | tee -a "$LOGFILE"
        return 1
    fi

    # Step 3: Compute bias + p-values
    echo "[JOB ${JOB_ID}] Step 3: Computing bias + p-values..." | tee -a "$LOGFILE"
    python "${BASE}/scripts/script_bias_cal.py" \
        --SpecMat "$SPECMAT" \
        --gw "$GW" \
        --Bias_Null "$NULL_BIAS" \
        --Bias_Out "$BIAS_OUT" \
        >> "$LOGFILE" 2>&1

    if [ $? -ne 0 ]; then
        echo "[JOB ${JOB_ID}] FAILED at Step 3" | tee -a "$LOGFILE"
        return 1
    fi

    echo "[JOB ${JOB_ID}] COMPLETED: ${TOPDIR}/${ATYPE}" | tee -a "$LOGFILE"
    echo "[JOB ${JOB_ID}] $(date)" | tee -a "$LOGFILE"
    return 0
}

export -f run_pipeline
export BASE GW SPECMAT_CENTER SPECMAT_QN LOGDIR

echo "============================================"
echo "Rerunning DDD_61 pipelines (32 jobs total)"
echo "Parallel: 3 at a time, 10 processes each"
echo "Logs: ${LOGDIR}"
echo "Started: $(date)"
echo "============================================"

# Define all jobs as tab-separated: TOPDIR ATYPE SPECMAT STEP1_ARGS
# We'll use GNU parallel or a simple batch approach

# Job definitions (one per line)
JOBS=(
    # --- Batch 1: Random and simple modes ---
    "systematic_random|Centering|${SPECMAT_CENTER}|--sampling_mode random --seed 42"
    "random|Centering|${SPECMAT_CENTER}|--sampling_mode random --seed 42"
    "random|QuantileNorm|${SPECMAT_QN}|--sampling_mode random --seed 42"

    # --- Batch 2: Best-of-N with standard vars ---
    "systematic_best_of_n_WB_mean_phastCons_n_CDS_bases_Best1000|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables n_CDS_bases,WB,mean_phastCons --use_best_of_n --n_candidates 1000 --use_percentile --seed 42"
    "set_level_matched|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables n_CDS_bases,WB,mean_phastCons --use_best_of_n --n_candidates 1000 --use_percentile --seed 42"
    "set_level_matched|QuantileNorm|${SPECMAT_QN}|--sampling_mode set_level_matched --matching_variables n_CDS_bases,WB,mean_phastCons --use_best_of_n --n_candidates 1000 --use_percentile --seed 42"

    # --- Batch 3: Best-of-N variants ---
    "conservation_model_WB_mean_phastCons_n_CDS_bases_Best1000|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables WB,mean_phastCons,n_CDS_bases --use_best_of_n --n_candidates 1000 --use_percentile --seed 42"
    "conservation_model_LOEUF_WB_n_CDS_bases_Best1000|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables LOEUF,WB,n_CDS_bases --use_best_of_n --n_candidates 1000 --use_percentile --seed 42"
    "set_level_matched_CDS_WB_LOEUF_Best1000|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables n_CDS_bases,WB,LOEUF --use_best_of_n --n_candidates 1000 --use_percentile --seed 42"

    # --- Batch 4: WB+CDS best-of-N variants ---
    "set_level_matched_WB_CDS|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables WB,n_CDS_bases --use_best_of_n --n_candidates 1000 --use_percentile --seed 42"
    "set_level_matched_best1000_WB_CDS|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables WB,n_CDS_bases --use_best_of_n --n_candidates 1000 --use_percentile --seed 42"
    "set_level_matched_best500_WB_CDS|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables WB,n_CDS_bases --use_best_of_n --n_candidates 500 --use_percentile --seed 42"

    # --- Batch 5: CDS-only and gene-by-gene ---
    "set_level_matched_CDS|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables n_CDS_bases --use_best_of_n --n_candidates 1000 --use_percentile --seed 42"
    "matched_WB_CDS|Centering|${SPECMAT_CENTER}|--sampling_mode matched --matching_variables WB,n_CDS_bases --kernel tricubic --bandwidth 10.0 --seed 42"
    "systematic_gene_by_gene_WB_mean_phastCons_n_CDS_bases_Tricubic|Centering|${SPECMAT_CENTER}|--sampling_mode matched --matching_variables n_CDS_bases,WB,mean_phastCons --kernel tricubic --bandwidth 100.0 --seed 42"

    # --- Batch 6: Rejection and SIS ---
    "systematic_rejection_WB_mean_phastCons_n_CDS_bases_Rejection|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables n_CDS_bases,WB,mean_phastCons --max_distance 0.15 --max_attempts 1000 --seed 42"
    "systematic_sis_WB_mean_phastCons_n_CDS_bases_SIS|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables n_CDS_bases,WB,mean_phastCons --use_sis --use_percentile --temperature 1.0 --adaptive_temp --seed 42"
    "set_level_matched_sis|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables n_CDS_bases,WB,mean_phastCons --use_sis --use_percentile --temperature 1.0 --adaptive_temp --seed 42"

    # --- Batch 7: SIS WB_CDS and propensity standard ---
    "set_level_matched_sis_WB_CDS|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables WB,n_CDS_bases --use_sis --use_percentile --temperature 1.0 --adaptive_temp --seed 42"
    "systematic_propensity_WB_mean_phastCons_n_CDS_bases_PropWeight_Tricubic|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables n_CDS_bases,WB,mean_phastCons --use_propensity --propensity_kernel tricubic --propensity_bandwidth 60.0 --add_noise --noise_scale 10.0 --relaxed_matching --loeuf_weight 1 --use_percentile --seed 42"
    "conservation_model_WB_mean_phastCons_n_CDS_bases_PropWeight_Tricubic|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables WB,mean_phastCons,n_CDS_bases --use_propensity --propensity_kernel tricubic --propensity_bandwidth 60.0 --add_noise --noise_scale 10.0 --relaxed_matching --loeuf_weight 1 --use_percentile --seed 42"

    # --- Batch 8: Propensity variants with LOEUF ---
    "constraint_model_LOEUF_WB_n_CDS_bases_PropWeight_Tricubic|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables LOEUF,WB,n_CDS_bases --use_propensity --propensity_kernel tricubic --propensity_bandwidth 60.0 --add_noise --noise_scale 10.0 --relaxed_matching --loeuf_weight 1 --use_percentile --seed 42"
    "set_level_matched_CDS_WB_LOEUF_PropWeight_Tricubic|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables n_CDS_bases,WB,LOEUF --use_propensity --propensity_kernel tricubic --propensity_bandwidth 60.0 --add_noise --noise_scale 10.0 --relaxed_matching --loeuf_weight 1 --use_percentile --seed 42"
    "set_level_matched_WB_CDS_LOEUF_PropWeight_Tricubic|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables WB,n_CDS_bases,LOEUF --use_propensity --propensity_kernel tricubic --propensity_bandwidth 60.0 --add_noise --noise_scale 10.0 --relaxed_matching --loeuf_weight 1 --use_percentile --seed 42"

    # --- Batch 9: More propensity variants ---
    "set_level_matched_CDS_LOEUF_PropWeight_Tricubic|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables n_CDS_bases,LOEUF --use_propensity --propensity_kernel tricubic --propensity_bandwidth 60.0 --add_noise --noise_scale 10.0 --relaxed_matching --loeuf_weight 1 --use_percentile --seed 42"
    "set_level_matched_WB_LOEUF_PropWeight_Tricubic|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables WB,LOEUF --use_propensity --propensity_kernel tricubic --propensity_bandwidth 60.0 --add_noise --noise_scale 10.0 --relaxed_matching --loeuf_weight 1 --use_percentile --seed 42"
    "set_level_matched_WB_mean_phastCons_n_CDS_bases_PropWeight_Tricubic|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables WB,mean_phastCons,n_CDS_bases --use_propensity --propensity_kernel tricubic --propensity_bandwidth 60.0 --add_noise --noise_scale 10.0 --relaxed_matching --loeuf_weight 1 --use_percentile --seed 42"

    # --- Batch 10: WB+CDS propensity variants ---
    "set_level_matched_WB_CDS_PropWeight|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables WB,n_CDS_bases --use_propensity --propensity_bandwidth 60.0 --add_noise --noise_scale 10.0 --relaxed_matching --loeuf_weight 1 --use_percentile --seed 42"
    "set_level_matched_WB_CDS_PropWeight_Tricubic|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables WB,n_CDS_bases --use_propensity --propensity_kernel tricubic --propensity_bandwidth 60.0 --add_noise --noise_scale 10.0 --relaxed_matching --loeuf_weight 1 --use_percentile --seed 42"
    "set_level_matched_WB_CDS_PropWeight_Uniform|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables WB,n_CDS_bases --use_propensity --propensity_kernel uniform --propensity_bandwidth 60.0 --add_noise --noise_scale 10.0 --relaxed_matching --loeuf_weight 1 --use_percentile --seed 42"

    # --- Batch 11: CDS+WB propensity tricubic v1 and v2 ---
    "set_level_matched_CDS_WB_PropWeight_Tricubic|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables n_CDS_bases,WB --use_propensity --propensity_kernel tricubic --propensity_bandwidth 60.0 --add_noise --noise_scale 10.0 --relaxed_matching --loeuf_weight 1 --use_percentile --seed 42"
    "set_level_matched_CDS_WB_PropWeight_Tricubic_v2|Centering|${SPECMAT_CENTER}|--sampling_mode set_level_matched --matching_variables n_CDS_bases,WB --use_propensity --propensity_kernel tricubic --propensity_bandwidth 60.0 --add_noise --noise_scale 10.0 --relaxed_matching --loeuf_weight 1 --use_percentile --seed 42"
)

TOTAL=${#JOBS[@]}
echo "Total jobs to run: ${TOTAL}"
echo ""

# Run in batches of 3
BATCH_SIZE=3
JOB_IDX=0

for ((i=0; i<TOTAL; i+=BATCH_SIZE)); do
    BATCH_END=$((i + BATCH_SIZE))
    if [ $BATCH_END -gt $TOTAL ]; then
        BATCH_END=$TOTAL
    fi

    echo "============================================"
    echo "BATCH: Jobs $((i+1)) to ${BATCH_END} of ${TOTAL}"
    echo "$(date)"
    echo "============================================"

    PIDS=()
    for ((j=i; j<BATCH_END; j++)); do
        IFS='|' read -r TOPDIR ATYPE SPECMAT STEP1_ARGS <<< "${JOBS[$j]}"
        JOB_NUM=$((j+1))

        run_pipeline "$TOPDIR" "$ATYPE" "$SPECMAT" "$STEP1_ARGS" "$JOB_NUM" &
        PIDS+=($!)
        echo "  Launched JOB ${JOB_NUM}: ${TOPDIR}/${ATYPE} (PID ${PIDS[-1]})"
    done

    # Wait for all jobs in this batch
    for PID in "${PIDS[@]}"; do
        wait $PID || {
            echo "WARNING: Process $PID failed"
            FAILED_JOBS="${FAILED_JOBS} PID:${PID}"
        }
    done

    echo "Batch complete: $(date)"
    echo ""
done

echo "============================================"
echo "ALL DONE: $(date)"
if [ -n "$FAILED_JOBS" ]; then
    echo "WARNING: Some jobs failed: ${FAILED_JOBS}"
    echo "Check logs in ${LOGDIR}"
else
    echo "All 32 jobs completed successfully"
fi
echo "============================================"
