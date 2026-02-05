#SPECMAT="/home/jw3514/Work/CellType_Psy/dat/Test.BiasMat/HumanCT.Spec.clip.noLowExp.cut1e4.csv"
#SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.csv"
#SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.0.Filt.Spec.clip.lowexp.cut1e4.csv"
SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.csv"
# Base output directory for bias results
RES_DIR="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/Spec_Bias_Jul07/"
OUT_BASE="$RES_DIR/CTRL/RandGeneBias_Spec"
mkdir -p "$RES_DIR"
mkdir -p "$OUT_BASE"

# Map keys to GW_DIR and output subdir suffix
declare -A GW_DIR_MAP
declare -A OUT_SUFFIX_MAP
declare -A DESC_MAP

GW_DIR_MAP["22q_small_del"]="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/RandGeneWeights/RandGene_22q_small_del/"
OUT_SUFFIX_MAP["22q_small_del"]="RandGene_22q_small_del"
DESC_MAP["22q_small_del"]="22q11.2 small deletion genes"

GW_DIR_MAP["22q_del"]="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/RandGeneWeights/RandGene_22q_del/"
OUT_SUFFIX_MAP["22q_del"]="RandGene_22q_del"
DESC_MAP["22q_del"]="22q11.2 deletion genes"

GW_DIR_MAP["ASD_HIQ"]="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/RandGeneWeights/RandGene_ASD_HIQ.top61/"
OUT_SUFFIX_MAP["ASD_HIQ"]="RandGene_ASD_HIQ.top61"
DESC_MAP["ASD_HIQ"]="HIQ Top60 ASD genes"

GW_DIR_MAP["ASD_LIQ"]="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/RandGeneWeights/RandGene_ASD_LIQ.top61/"
OUT_SUFFIX_MAP["ASD_LIQ"]="RandGene_ASD_LIQ.top61"
DESC_MAP["ASD_LIQ"]="LIQ Top60 ASD genes"

GW_DIR_MAP["SCZ"]="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/RandGeneWeights/RandGene_SCZ.top61/"
OUT_SUFFIX_MAP["SCZ"]="RandGene_SCZ.top61"
DESC_MAP["SCZ"]="SCZ Top60 genes"

GW_DIR_MAP["DDD"]="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/RandGeneWeights/RandGene_DDD.top61/"
OUT_SUFFIX_MAP["DDD"]="RandGene_DDD.top61"
DESC_MAP["DDD"]="DDD Top60 genes"

GW_DIR_MAP["UKBB_VNR"]="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/RandGeneWeights/RandGene_UKBB_VNR.top61/"
OUT_SUFFIX_MAP["UKBB_VNR"]="RandGene_UKBB_VNR.top61"
DESC_MAP["UKBB_VNR"]="UKBB VNR genes"

GW_DIR_MAP["UKBB_EDU"]="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/RandGeneWeights/RandGene_UKBB_EDU.top61/"
OUT_SUFFIX_MAP["UKBB_EDU"]="RandGene_UKBB_EDU.top61"
DESC_MAP["UKBB_EDU"]="UKBB EDU genes"

GW_DIR_MAP["UKBB_RT"]="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/RandGeneWeights/RandGene_UKBB_RT.top61/"
OUT_SUFFIX_MAP["UKBB_RT"]="RandGene_UKBB_RT.top61"
DESC_MAP["UKBB_RT"]="UKBB RT genes"

for key in "22q_small_del" "22q_del" "ASD_HIQ" "ASD_LIQ" "SCZ" "DDD" "UKBB_VNR" "UKBB_EDU" "UKBB_RT"; do
    echo "Running ${DESC_MAP[$key]}"
    GW_DIR="${GW_DIR_MAP[$key]}"
    BiasOUTDIR="$OUT_BASE/${OUT_SUFFIX_MAP[$key]}/"
    mkdir -p "$BiasOUTDIR"
    python script_run_ctrl_sim.v3.py -m human_ct_bias -o "$BiasOUTDIR" --SpecMat "$SPECMAT" --GW_Dir "$GW_DIR"
done

