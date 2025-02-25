# Null distribution under different null models for ASD, HIQ ASD and LIQ ASD 

#BiasMat="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip5.Z2.clip3.Dec30.csv"
BiasMat="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip5.Z2.clip3.Jan21.csv"



###########
# ASD
###########
# ASD_LIQ_MUT_FILE="/home/jw3514/Work/CellType_Psy/dat/GeneWeights/ASD.LowIQ.Mutable.csv"
# ASD_HIQ_MUT_FILE="/home/jw3514/Work/CellType_Psy/dat/GeneWeights/ASD.HighIQ.Mutable.csv"
# ASD_ALL_MUT_FILE="/home/jw3514/Work/CellType_Psy/dat/GeneWeights/ASD.All.Mutable.csv"
# ASD_ALL_GENE_FILE="/home/jw3514/Work/CellType_Psy/dat/GeneWeights/ASD.ExomeWide.Denovo.csv"

# #1. Random mutation 
# OutDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/WeightRedistribute_Jan17/ASD_HIQ_Type1/"
# python script_WeightRedistribute.py --dataset ASD --sim_type 1 --mutable $ASD_HIQ_MUT_FILE --BiasMat $BiasMat -d $OutDIR --n_processes 30

# OutDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/WeightRedistribute_Jan17/ASD_LIQ_Type1/"
# python script_WeightRedistribute.py --dataset ASD --sim_type 1 --mutable $ASD_LIQ_MUT_FILE --BiasMat $BiasMat -d $OutDIR --n_processes 30

# OutDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/WeightRedistribute_Jan17/ASD_ALL_Type1/"
# python script_WeightRedistribute.py --dataset ASD --sim_type 1 --mutable $ASD_ALL_MUT_FILE --BiasMat $BiasMat -d $OutDIR --n_processes 30

# OutDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/WeightRedistribute_Jan17/ASD_ExomeWide_Type1/"
# python script_WeightRedistribute.py --dataset ASD_ALL --sim_type 1 --mutable $ASD_ALL_GENE_FILE --BiasMat $BiasMat -d $OutDIR --n_processes 30

# # 2. Random mutation according to background mutation rate
# OutDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/WeightRedistribute_Jan17/ASD_HIQ_Type2/"
# python script_WeightRedistribute.py --dataset ASD --sim_type 2 --mutable $ASD_HIQ_MUT_FILE --BiasMat $BiasMat -d $OutDIR --n_processes 30

# OutDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/WeightRedistribute_Jan17/ASD_LIQ_Type2/"
# python script_WeightRedistribute.py --dataset ASD --sim_type 2 --mutable $ASD_LIQ_MUT_FILE --BiasMat $BiasMat -d $OutDIR --n_processes 30

# OutDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/WeightRedistribute_Jan17/ASD_ALL_Type2/"
# python script_WeightRedistribute.py --dataset ASD --sim_type 2 --mutable $ASD_ALL_MUT_FILE --BiasMat $BiasMat -d $OutDIR --n_processes 30

# OutDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/WeightRedistribute_Jan17/ASD_ExomeWide_Type2/"
# python script_WeightRedistribute.py --dataset ASD_ALL --sim_type 2 --mutable $ASD_ALL_GENE_FILE --BiasMat $BiasMat -d $OutDIR --n_processes 30

# # # 3. Random mutation according to background mutation rate and assign new weights to random genes
# OutDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/WeightRedistribute_Jan17/ASD_HIQ_Type3/"
# python script_WeightRedistribute.py --dataset ASD --sim_type 3 --mutable $ASD_HIQ_MUT_FILE --BiasMat $BiasMat -d $OutDIR --n_processes 30

# OutDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/WeightRedistribute_Jan17/ASD_LIQ_Type3/"
# python script_WeightRedistribute.py --dataset ASD --sim_type 3 --mutable $ASD_LIQ_MUT_FILE --BiasMat $BiasMat -d $OutDIR --n_processes 30

# OutDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/WeightRedistribute_Jan17/ASD_ALL_Type3/"
# python script_WeightRedistribute.py --dataset ASD --sim_type 3 --mutable $ASD_ALL_MUT_FILE --BiasMat $BiasMat -d $OutDIR --n_processes 30

# OutDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/WeightRedistribute_Jan17/ASD_ExomeWide_Type3/"
# python script_WeightRedistribute.py --dataset ASD_ALL --sim_type 3 --mutable $ASD_ALL_GENE_FILE --BiasMat $BiasMat -d $OutDIR --n_processes 30

###########
# SCZ
###########

#1. Random mutation 
SCZ_MUT_FILE="/home/jw3514/Work/CellType_Psy/dat/GeneWeights/SCZ.top61.Mutable.csv"
OutDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/WeightRedistribute_Jan17/SCZ_Type1/"
python script_WeightRedistribute.py --dataset SCZ --sim_type 1 --mutable $SCZ_MUT_FILE --BiasMat $BiasMat -d $OutDIR --n_processes 30

# # 2. Random mutation according to control mutation counts of each gene 
# SCZ_MUT_FILE="/home/jw3514/Work/CellType_Psy/dat/GeneWeights/SCZ.top61.Mutable.csv"
# OutDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/WeightRedistribute_Jan17/SCZ_Type2/"
# python script_WeightRedistribute.py --dataset SCZ --sim_type 2 --mutable $SCZ_MUT_FILE --BiasMat $BiasMat -d $OutDIR --n_processes 30

# # 3. Random mutation according to control mutation counts of each gene 
# SCZ_MUT_FILE="/home/jw3514/Work/CellType_Psy/dat/GeneWeights/SCZ.top61.Mutable.csv"
# OutDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/WeightRedistribute_Jan17/SCZ_Type3/"
# python script_WeightRedistribute.py --dataset SCZ --sim_type 3 --mutable $SCZ_MUT_FILE --BiasMat $BiasMat -d $OutDIR --n_processes 30

# # 4. Random mutation according to control mutation counts of each gene 
# SCZ_MUT_FILE="/home/jw3514/Work/CellType_Psy/dat/GeneWeights/SCZ.top61.Mutable.csv"
# OutDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/WeightRedistribute_Jan17/SCZ_Type4/"
# python script_WeightRedistribute.py --dataset SCZ --sim_type 4 --mutable $SCZ_MUT_FILE --BiasMat $BiasMat -d $OutDIR --n_processes 30

###########
# NDD
###########
# NDD_MUT_FILE="/home/jw3514/Work/CellType_Psy/dat/GeneWeights/DDD.top61.GeneList.csv"
# #1. Random mutation 
# OutDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/WeightRedistribute_Jan17/NDD_Type1/"
# python script_WeightRedistribute.py --dataset NDD --sim_type 1 --mutable $NDD_MUT_FILE --BiasMat $BiasMat -d $OutDIR --n_processes 30

# # 2. Random mutation according to background mutation rate
# OutDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/WeightRedistribute_Jan17/NDD_Type2/"
# python script_WeightRedistribute.py --dataset NDD --sim_type 2 --mutable $NDD_MUT_FILE --BiasMat $BiasMat -d $OutDIR --n_processes 30

# # 3. Random mutation according to background mutation rate and assign new weights to random genes
# OutDIR="/home/jw3514/Work/CellType_Psy/dat/CTRL/WeightRedistribute_Jan17/NDD_Type3/"
# python script_WeightRedistribute.py --dataset NDD --sim_type 3 --mutable $NDD_MUT_FILE --BiasMat $BiasMat -d $OutDIR --n_processes 30