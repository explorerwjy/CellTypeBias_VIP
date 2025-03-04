SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip5.Z2.clip3.Dec30.csv"

###########################
# 22q11.2 genes 
###########################

# simulate gene weights
# GW4SIM="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/GeneWeights/X22q.gw.csv"
# GW_DIR="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/RandGeneWeights/RandGenes_22q/"
# N_SIMS=10000
# mkdir -p $GW_DIR
# python script_run_ctrl_sim.v2.py -m gw -o $GW_DIR --n_sims $N_SIMS  -w $GW4SIM
# # gene permutation test
# BiasOUTDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/RandGeneBias/RandGenes_22q/"
# mkdir -p $BiasOUTDIR
# python script_run_ctrl_sim.v2.py -m human_ct_bias -o $BiasOUTDIR --SpecMat $SPECMAT --GW_Dir $GW_DIR

###########################
# 22q11 small deletion genes 
###########################
# simulate gene weights
# GW4SIM="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/GeneWeights/X22q.mousemodel.gw.csv"
# GW_DIR="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/RandGeneWeights/RandGenes_22q_small_del/"
# N_SIMS=10000
# mkdir -p $GW_DIR
# python script_run_ctrl_sim.v2.py -m gw -o $GW_DIR --n_sims $N_SIMS  -w $GW4SIM
# # gene permutation test
# BiasOUTDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/RandGeneBias/RandGenes_22q_small_del/"
# mkdir -p $BiasOUTDIR
# python script_run_ctrl_sim.v2.py -m human_ct_bias -o $BiasOUTDIR --SpecMat $SPECMAT --GW_Dir $GW_DIR



###########################
# LIQ Top60 (51) ASD genes 
###########################
# simulate gene weights
GW4SIM="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/GeneWeights/HIQ.top61.nopLI.gw"
GW_DIR="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/RandGeneWeights/RandGenes_ASD_HIQ.top61/"
N_SIMS=10000
mkdir -p $GW_DIR
python script_run_ctrl_sim.v2.py -m gw -o $GW_DIR --n_sims $N_SIMS  -w $GW4SIM

BiasOUTDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/RandGeneBias/RandGenes_ASD_HIQ.top61/"
mkdir -p $BiasOUTDIR
python script_run_ctrl_sim.v2.py -m human_ct_bias -o $BiasOUTDIR --SpecMat $SPECMAT --GW_Dir $GW_DIR

###########################
# HIQ Top60 (51) ASD genes 
###########################
# # simulate gene weights
GW4SIM="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/GeneWeights/LIQ.top61.nopLI.gw"
GW_DIR="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/RandGeneWeights/RandGenes_ASD_LIQ.top61/"
N_SIMS=10000
mkdir -p $GW_DIR
python script_run_ctrl_sim.v2.py -m gw -o $GW_DIR --n_sims $N_SIMS  -w $GW4SIM

BiasOUTDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/RandGeneBias/RandGenes_ASD_LIQ.top61/"
mkdir -p $BiasOUTDIR
python script_run_ctrl_sim.v2.py -m human_ct_bias -o $BiasOUTDIR --SpecMat $SPECMAT --GW_Dir $GW_DIR

###########################
# DD/ID genes 
###########################
# # simulate gene weights
# GW4SIM="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/GeneWeights/DDD.top61.gw.csv"
# GW_DIR="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/RandGeneWeights/RandGenes_DDD.top61/"
# N_SIMS=10000
# mkdir -p $GW_DIR
# python script_run_ctrl_sim.v2.py -m gw -o $GW_DIR --n_sims $N_SIMS  -w $GW4SIM

# BiasOUTDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/RandGeneBias/RandGenes_DDD.top61/"
# mkdir -p $BiasOUTDIR
# python script_run_ctrl_sim.v2.py -m human_ct_bias -o $BiasOUTDIR --SpecMat $SPECMAT --GW_Dir $GW_DIR

###########################
# VNR genes 
###########################
# # simulate gene weights
# GW4SIM="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/GeneWeights/UKBB_VNR_Neg_GW_61.csv"
# GW_DIR="/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/RandGeneWeights/RandGenes_UKBB_VNR.top61/"
# N_SIMS=10000
# mkdir -p $GW_DIR
# python script_run_ctrl_sim.v2.py -m gw -o $GW_DIR --n_sims $N_SIMS  -w $GW4SIM

# BiasOUTDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/RandGeneBias/RandGenes_UKBB_VNR.top61/"
# mkdir -p $BiasOUTDIR
# python script_run_ctrl_sim.v2.py -m human_ct_bias -o $BiasOUTDIR --SpecMat $SPECMAT --GW_Dir $GW_DIR

#####################################################
# Other gene set not included in the paper
#####################################################

###########################
# LIQ Top500 (163) ASD genes 
###########################
# simulate gene weights
# GW4SIM="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/Rare/LIQ.top163.gw"
# GW_DIR="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/CTRL/RandG_w_Weights/RandGenes_ASD_LIQ_top163/"
# N_SIMS=10000
# mkdir -p $GW_DIR
# python script_run_ctrl_sim.v2.py -m gw -o $GW_DIR --n_sims $N_SIMS  -w $GW4SIM

# BiasOUTDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/RandGenes/RandGene_ASD_LIQ_top163_clip3.3_Feb20/"
# SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip5.Z2.clip3.Jan21.csv"
# GW_DIR="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/CTRL/RandG_w_Weights/RandGenes_ASD_LIQ_top163/"
# mkdir -p $BiasOUTDIR
# python script_run_ctrl_sim.v2.py -m human_ct_bias -o $BiasOUTDIR --SpecMat $SPECMAT --GW_Dir $GW_DIR

###########################
# HIQ Top500 (163) ASD genes 
###########################
# simulate gene weights
# GW4SIM="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/Rare/HIQ.top188.gw"
# GW_DIR="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/CTRL/RandG_w_Weights/RandGenes_ASD_HIQ_top188/"
# N_SIMS=10000
# mkdir -p $GW_DIR
# python script_run_ctrl_sim.v2.py -m gw -o $GW_DIR --n_sims $N_SIMS  -w $GW4SIM

# BiasOUTDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/RandGenes/RandGene_ASD_HIQ_top188_clip3.3_Feb20/"
# SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip5.Z2.clip3.Jan21.csv"
# GW_DIR="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/CTRL/RandG_w_Weights/RandGenes_ASD_HIQ_top188/"
# mkdir -p $BiasOUTDIR
# python script_run_ctrl_sim.v2.py -m human_ct_bias -o $BiasOUTDIR --SpecMat $SPECMAT --GW_Dir $GW_DIR

###########################
# Dec30 Z2 Matrix Z1 clip 5, 60 genes weight transfer 
###########################
# Testing Dec30 Z2 Matrix Z1 clip 5, 60 genes weight transfer 
# SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip5.Z2.clip3.Dec30.csv"
# OUTDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/RandGenes/RandGene_60_Z1clip5.Z2.clip3.Dec30/"
# GW_DIR="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/CTRL/RandG_w_Weights/RandGenes_60/"
# mkdir -p $OUTDIR
#nohup python script_run_ctrl_sim.v2.py -m human_ct_bias -o $OUTDIR --SpecMat $SPECMAT --GW_Dir $GW_DIR > log.txt &

###########################
###########################
# GW4SIM="/home/jw3514/Work/CellType_Psy/dat3/SCZ_MutCount_61.gw"
# OUTDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/ExpMatch_Ctrl/SCZ_MutCount_61/"
# MatchDir="/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/MatchGenes.HumanCT/"
# N_SIMS=10000
# mkdir -p $OUTDIR
#python script_ctrl_sim_exp_match.py -m gw -o $OUTDIR --n_sims $N_SIMS  -w $GW4SIM --MatchDir $MatchDir

# SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip3.Z2.Jan08.Z2clip3.csv"
# OUTDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/ExpMatchGenes/SCZ_MutCount_61/"
# GW_DIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/ExpMatch_Ctrl/SCZ_MutCount_61/weighted/"
# mkdir -p $OUTDIR
#nohup python script_ctrl_sim_exp_match.py -m human_ct_bias -o $OUTDIR --SpecMat $SPECMAT --GW_Dir $GW_DIR > log.txt &

# Testing Dec30 Z2 Matrix Z1 clip 5, 60 genes weight transfer 
# SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip5.Z2.clip3.Dec30.csv"
# OUTDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/RandGenes/RandGene_60_Z1clip5.Z2.clip3.Dec30/"
# GW_DIR="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/CTRL/RandG_w_Weights/RandGenes_60/"
# mkdir -p $OUTDIR
#nohup python script_run_ctrl_sim.v2.py -m human_ct_bias -o $OUTDIR --SpecMat $SPECMAT --GW_Dir $GW_DIR > log.txt &

# Testing Dec30 Z2 Matrix Z2 clip at 3
# SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip5.Z2.clip3.Dec30.csv"
# OUTDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/RandGenes/UniformWeights/RandGene_60_Z1clip5.Z2.clip3.Dec30/"
# GW_DIR="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/CTRL/RandG_UniformWeights/RandGenes_60/"
# mkdir -p $OUTDIR
#nohup python script_run_ctrl_sim.v2.py -m human_ct_bias -o $OUTDIR --SpecMat $SPECMAT --GW_Dir $GW_DIR > log.txt &

# Testing Dec30 Z2 Matrix Z1 clip 5, 60 genes weight transfer 
# SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip5.Z2.Dec30.csv"
# OUTDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/RandGenes/RandGene_60_Z1clip5.Z2.Dec30/"
# GW_DIR="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/CTRL/RandG_w_Weights/RandGenes_60/"
# mkdir -p $OUTDIR
#nohup python script_run_ctrl_sim.v2.py -m human_ct_bias -o $OUTDIR --SpecMat $SPECMAT --GW_Dir $GW_DIR > log.txt &

# # Testing Dec30 Z2 Matrix lowexp_1e-3
# SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip5.lowexp_1e-3.Z2.Dec30.csv"
# OUTDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/RandGenes/UniformWeights/RandGene_60_Z1clip5.lowexp_1e-3.Z2.Dec30/"
# #GW_DIR="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/CTRL/RandG_w_Weights/RandGenes_60/"
# GW_DIR="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/CTRL/RandG_UniformWeights/RandGenes_60/"
# mkdir -p $OUTDIR
# #nohup python script_run_ctrl_sim.v2.py -m human_ct_bias -o $OUTDIR --SpecMat $SPECMAT --GW_Dir $GW_DIR > log.txt &

# # Testing Dec30 Z2 Matrix 
# SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip5.Z2.Dec30.csv"
# OUTDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/RandGenes/UniformWeights/RandGene_60_Z1clip5.Z2.Dec30/"
# #GW_DIR="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/CTRL/RandG_w_Weights/RandGenes_60/"
# GW_DIR="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/CTRL/RandG_UniformWeights/RandGenes_60/"
# mkdir -p $OUTDIR
# #nohup python script_run_ctrl_sim.v2.py -m human_ct_bias -o $OUTDIR --SpecMat $SPECMAT --GW_Dir $GW_DIR > log.txt &

# SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip5.Z2.Dec30.csv"
# GW_DIR="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/CTRL/RandG_w_Weights/RandGenes_46_22q"
# OUTDIR="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/CTRL/RandG_w_Weights/RandGenes_46_22q_Dec30/"
# mkdir -p $OUTDIR
# #nohup python script_run_ctrl_sim.v2.py -m human_ct_bias -o $OUTDIR --SpecMat $SPECMAT --GW_Dir $GW_DIR > log.txt &


# # Testing Nov14 Z2 Matrix 
# SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.log2.Z2.HCT.z1clip3.csv"
# OUTDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/RandGenes/UniformWeights/RandGene_60_z1clip3.3"
# #GW_DIR="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/CTRL/RandG_w_Weights/RandGenes_60/"
# GW_DIR="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/CTRL/RandG_UniformWeights/RandGenes_60/"
# mkdir -p $OUTDIR
# #python script_run_ctrl_sim.v2.py -m human_ct_bias -o $OUTDIR --SpecMat $SPECMAT --GW_Dir $GW_DIR

# # Random Genes 

# SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.log2.Z2.HCT.z1clip3.csv"
# OUTDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/RandGenes/UniformWeights/RandGene_60_z1clip3.3"
# #GW_DIR="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/CTRL/RandG_w_Weights/RandGenes_60/"
# GW_DIR="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/CTRL/RandG_UniformWeights/RandGenes_60/"
# mkdir -p $OUTDIR
# #python script_run_ctrl_sim.v2.py -m human_ct_bias -o $OUTDIR --SpecMat $SPECMAT --GW_Dir $GW_DIR

# #SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Subcluster.log2.Z2.csv"
# SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Residue.Z2.csv"
# GW4SIM="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/Rare/GW.ASD.top60.txt"
# OUTDIR="/home/jw3514/Work/UNIMED/dat/Genetics/GeneWeights/CTRL/RandG_w_Weights/RandGenes_60_ResidueZ2/"
# N_SIMS=10000
# mkdir -p $OUTDIR
# #python script_run_ctrl_sim.v2.py -m gw -o $OUTDIR --n_sims $N_SIMS  -w $GW4SIM --SpecMat $SPECMAT

# BiasOUTDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/RandGenes/RandGenes_60_ResidueZ2/"
# SPECMAT="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Residue.Z2.csv"
# GW_DIR=$OUTDIR
# mkdir -p $BiasOUTDIR
# #python script_run_ctrl_sim.v2.py -m human_ct_bias -o $BiasOUTDIR --SpecMat $SPECMAT --GW_Dir $GW_DIR




