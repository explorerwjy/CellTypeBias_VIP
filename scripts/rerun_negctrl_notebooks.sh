#!/bin/bash
# Re-execute the 6 notebooks that consume NegCtrl bias_addP files,
# now that production NegCtrls have been switched from uniform to abs(β).
# Per-notebook logging to logs/nbrun_<name>.log; aggregate status to stdout.

set -u
cd /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP
mkdir -p logs

NOTEBOOKS=(
    "notebooks_rebuttal/NegativeControl_BiasPlot"
    "notebooks_rebuttal/NegativeControl_Investigation"
    "notebooks_rebuttal/QQ_Plot_Systematic"
    "notebooks/Figures_Supp"
    "notebooks/Bias_Mutation_Weights"
    "notebooks/GWAS_RareMutation_Comparison"
)

PASS=()
FAIL=()
for nb in "${NOTEBOOKS[@]}"; do
    name=$(basename "$nb")
    log="logs/nbrun_${name}.log"
    echo "[$(date +%H:%M:%S)] >>> $nb"
    {
        conda run -n gencic --no-capture-output jupytext --sync "${nb}.py" 2>&1
        conda run -n gencic --no-capture-output jupyter nbconvert --to notebook --execute --inplace "${nb}.ipynb" 2>&1
    } > "$log" 2>&1
    rc=$?
    if [ $rc -eq 0 ]; then
        echo "[$(date +%H:%M:%S)] OK    $nb"
        PASS+=("$nb")
    else
        echo "[$(date +%H:%M:%S)] FAIL  $nb (rc=$rc, see $log)"
        FAIL+=("$nb")
    fi
done

echo
echo "==================== SUMMARY ===================="
echo "PASSED (${#PASS[@]}): ${PASS[@]}"
echo "FAILED (${#FAIL[@]}): ${FAIL[@]}"
echo "================================================="
exit ${#FAIL[@]}
