SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.csv"
#SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.0.Filt.Spec.clip.lowexp.cut1e4.csv"
# Base output directory for bias results


# Base directories for gene weights and output
GW_BASE="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/GeneWeights"

RES_DIR="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/Spec_Bias_Jul07/"
OUT_BASE="$RES_DIR/CTRL/RandGeneWeights"
mkdir -p "$RES_DIR"
mkdir -p "$OUT_BASE"
N_SIMS=10000

###########################
# Simulate gene weights
###########################

declare -A GW4SIM_MAP
declare -A GW_DIR_MAP
declare -A DESC_MAP

GW4SIM_MAP["22q_small_del"]="$GW_BASE/X22q.mousemodel.gw.csv"
GW_DIR_MAP["22q_small_del"]="$OUT_BASE/RandGene_22q_small_del/"
DESC_MAP["22q_small_del"]="22q11.2 small deletion genes"

GW4SIM_MAP["22q_del"]="$GW_BASE/X22q.gw.csv"
GW_DIR_MAP["22q_del"]="$OUT_BASE/RandGene_22q_del/"
DESC_MAP["22q_del"]="22q11.2 deletion genes"

GW4SIM_MAP["ASD_HIQ"]="$GW_BASE/HIQ.top61.nopLI.LGD_Dmis_SameWeight.gw"
GW_DIR_MAP["ASD_HIQ"]="$OUT_BASE/RandGene_ASD_HIQ.top61/"
DESC_MAP["ASD_HIQ"]="HIQ ASD genes"

GW4SIM_MAP["ASD_LIQ"]="$GW_BASE/LIQ.top61.nopLI.LGD_Dmis_SameWeight.gw"
GW_DIR_MAP["ASD_LIQ"]="$OUT_BASE/RandGene_ASD_LIQ.top61/"
DESC_MAP["ASD_LIQ"]="LIQ ASD genes"

GW4SIM_MAP["SCZ"]="$GW_BASE/SCZ.top61.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw"
GW_DIR_MAP["SCZ"]="$OUT_BASE/RandGene_SCZ.top61/"
DESC_MAP["SCZ"]="SCZ genes"

GW4SIM_MAP["DDD"]="$GW_BASE/DDD.top61.gw.csv"
GW_DIR_MAP["DDD"]="$OUT_BASE/RandGene_DDD.top61/"
DESC_MAP["DDD"]="DDD genes"

GW4SIM_MAP["UKBB_VNR"]="$GW_BASE/UKBB_VNR_Neg_GW_61.csv"
GW_DIR_MAP["UKBB_VNR"]="$OUT_BASE/RandGene_UKBB_VNR.top61/"
DESC_MAP["UKBB_VNR"]="UKBB VNR genes"

GW4SIM_MAP["UKBB_EDU"]="$GW_BASE/UKBB_EDU_Neg_GW_61.csv"
GW_DIR_MAP["UKBB_EDU"]="$OUT_BASE/RandGene_UKBB_EDU.top61/"
DESC_MAP["UKBB_EDU"]="UKBB EDU genes"

GW4SIM_MAP["UKBB_RT"]="$GW_BASE/UKBB_RT_Neg_GW_61.csv"
GW_DIR_MAP["UKBB_RT"]="$OUT_BASE/RandGene_UKBB_RT.top61/"
DESC_MAP["UKBB_RT"]="UKBB RT genes"

for key in "22q_small_del" "22q_del" "ASD_HIQ" "ASD_LIQ" "SCZ" "DDD" "UKBB_VNR" "UKBB_EDU" "UKBB_RT"; do
    echo "Simulating gene weights for ${DESC_MAP[$key]}"
    mkdir -p "${GW_DIR_MAP[$key]}"
    python script_run_ctrl_sim.v3.py -m gw \
        -o "${GW_DIR_MAP[$key]}" \
        --n_sims $N_SIMS \
        -w "${GW4SIM_MAP[$key]}" \
        --SpecMat $SPECMAT
done
