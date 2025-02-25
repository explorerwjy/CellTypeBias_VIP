#!/bin/bash
Start=0
Step=500
End=20000
Script="script_Z2_calculation.py"

run_z2_calculation() {
    local ExpMat="$1"
    local MatchDir="$2"
    local OutDir="$3"
    local OutFil="$4"

    mkdir -p "$OutDir"
    parallel -j 20 bash run_Z2.sh {} $Step "$ExpMat" "$MatchDir" "$OutDir" ::: $(seq $Start $Step $End)
    python script_CombineZ2.py "$OutDir" "$OutFil"
}
# Example usage:

ExpMat="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Exp.Z1.clip3.csv"
MatchDir="/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/MatchGenes.HumanCT/"
OutDir="/home/jw3514/Work/CellType_Psy/dat/Z2.Split/Human.Cluster.Log2Mean.Z1clip3.Z2.Jan08/"
OutFil="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip3.Z2.Jan08.csv"
run_z2_calculation "$ExpMat" "$MatchDir" "$OutDir" "OutFil"
python script_CombineZ2.py "$OutDir" "$OutFil"

ExpMat="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Exp.Z1.clip3.csv"
MatchDir="/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/MatchGenes.HumanCT/"
OutDir="/home/jw3514/Work/CellType_Psy/dat/Z2.Split/Human.Cluster.Log2Mean.Z1clip3.Z2.Jan08/"
OutFil="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip3.Z2.Jan08.csv"
run_z2_calculation "$ExpMat" "$MatchDir" "$OutDir" "OutFil"
python script_CombineZ2.py "$OutDir" "$OutFil"




# ExpMat="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Spec.Test.Nov14.csv"
# MatchDir="/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/MatchGenes.HumanCT/"
# OutDir="/home/jw3514/Work/CellType_Psy/dat/Z2.Split/Human.Cluster.Spec.Test.Nov14.Z2/"
# OutFil="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Spec.Test.Nov14.Z2.csv"
# #run_z2_calculation "$ExpMat" "$MatchDir" "$OutDir" "OutFil"
# python script_CombineZ2.py "$OutDir" "$OutFil"

# Human Cluster level with Z1 clip 5 with low exp of 1e-3 (expression lower then 1e-3 give -5 z1), and normal Z2 . 
# ExpMat="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Exp.Z1.clip5.lowexp_1e-3.csv"
# MatchDir="/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/MatchGenes.HumanCT/"
# OutDir="/home/jw3514/Work/CellType_Psy/dat/Z2.Split/Human.Cluster.Log2Mean.Z1clip5.lowexp_1e-3.Z2.Dec30/"
# OutFil="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip5.lowexp_1e-3.Z2.Dec30.csv"
# run_z2_calculation "$ExpMat" "$MatchDir" "$OutDir" "OutFil"
# python script_CombineZ2.py "$OutDir" "$OutFil"

# Human Cluster level with Z1 clip 5 (expression 0 give -3 z1), and normal Z2 . 
# ExpMat="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Exp.Z1.clip3.csv"
# MatchDir="/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/MatchGenes.HumanCT/"
# OutDir="/home/jw3514/Work/CellType_Psy/dat/Z2.Split/Human.Cluster.Log2Mean.Z1clip3.Z2.Jan08/"
# OutFil="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip3.Z2.Jan08.csv"
# run_z2_calculation "$ExpMat" "$MatchDir" "$OutDir" "OutFil"
# python script_CombineZ2.py "$OutDir" "$OutFil"

# Human Cluster level with Z1 clip 5 (expression 0 give -3 z1), and normal Z2 . 
# ExpMat="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Exp.Z1.clip5.csv"
# MatchDir="/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/MatchGenes.HumanCT/"
# OutDir="/home/jw3514/Work/CellType_Psy/dat/Z2.Split/Human.Cluster.Log2Mean.Z1clip5.Z2.Dec30/"
# OutFil="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip5.Z2.Dec30.csv"
# run_z2_calculation "$ExpMat" "$MatchDir" "$OutDir" "OutFil"
# python script_CombineZ2.py "$OutDir" "$OutFil"

# Example usage:
# ExpMat="/home/jw3514/Work/CellType_Psy/dat/HumanCT.Residue.csv"
# MatchDir="/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/MatchGenes.HumanCT/"
# OutDir="/home/jw3514/Work/CellType_Psy/dat/Z2.Split/MouseSTR_ISH_ResidueZ2.Z2/"
# OutFil="/home/jw3514/Work/CellType_Psy/dat/ResidueBiasMat/MouseSTR_ISH_ResidueZ2.csv"
# run_z2_calculation "$ExpMat" "$MatchDir" "$OutDir"

#ExpMat="/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Subcluster.log2.Z1.clip3.csv"
#MatchDir="/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/MatchGenes.HumanCT/"
#OutDir="/home/jw3514/Work/CellType_Psy/dat2/Z2.Split/Z2.Subcluster.HumanCTMatch.clip3/"
#mkdir -p $OutDir
#parallel -j 20 bash run_Z2.sh {} $Step $ExpMat $MatchDir $OutDir ::: $(seq $Start $Step $End)

#ExpMat="/home/jw3514/Work/CellType_Psy/dat/Human.CT.Exp.Entrez.log2.Z1.clip.csv"
#MatchDir="/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/MatchGenes.HumanCT/"
#OutDir="/home/jw3514/Work/CellType_Psy/dat2/Z2.Split/Z2.ExpL.HumanCTMatch/"
#mkdir -p $OutDir
#parallel -j 20 bash run_Z2.sh {} $Step $ExpMat $MatchDir $OutDir ::: $(seq $Start $Step $End)

#ExpMat="/home/jw3514/Work/CellType_Psy/dat/Human.CT.Exp.Entrez.log2.Z1.clip.csv"
#MatchDir="/home/jw3514/Work/ASD_Circuits/dat/genes/BrainSpanMatch_uniform_kernal/"
#OutDir="/home/jw3514/Work/CellType_Psy/dat2/Z2.Split/Z2.ExpL.BrainSpanMatch/"
#mkdir -p $OutDir
#parallel -j 20 bash run_Z2.sh {} $Step $ExpMat $MatchDir $OutDir ::: $(seq $Start $Step $End)


#ExpMat="/home/jw3514/Work/CellType_Psy/dat/Human.CT.Exp.Entrez.log2.Z1.clip3.csv"
#MatchDir="/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/MatchGenes.HumanCT/"
#OutDir="/home/jw3514/Work/CellType_Psy/dat2/Z2.Split/Z2.ExpL.HumanCTMatch_z1clip3/"
#mkdir -p $OutDir
#parallel -j 20 bash run_Z2.sh {} $Step $ExpMat $MatchDir $OutDir ::: $(seq $Start $Step $End)

# ExpMat="/home/jw3514/Work/CellType_Psy/dat/Human.CT.Exp.Entrez.log2.Z1.clip.csv"
# MatchDir="/home/jw3514/Work/ASD_Circuits/dat/genes/BrainSpanMatch_uniform_kernal/"
# OutDir="/home/jw3514/Work/CellType_Psy/dat2/Z2.Split/Z2.ExpL.BrainSpanMatch/"
# mkdir -p $OutDir
#parallel -j 20 bash run_Z2.sh {} $Step $ExpMat $MatchDir $OutDir ::: $(seq $Start $Step $End)

# ExpMat="../dat/Human.CT.Exp.Entrez.csv"
# MatchDir="/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/MatchGenes.HumanCT/"
# OutDir="/home/jw3514/Work/CellType_Psy/dat2/Z2.Split/Zmatch.ExpL.HumanCTMatch/"
# mkdir -p $OutDir
#parallel -j 20 bash run_Z2.sh {} $Step $ExpMat $MatchDir $OutDir ::: $(seq $Start $Step $End)

# ExpMat="../dat/Human.CT.Exp.Entrez.csv"
# MatchDir="/home/jw3514/Work/ASD_Circuits/dat/genes/BrainSpanMatch_uniform_kernal/"
# OutDir="/home/jw3514/Work/CellType_Psy/dat2/Z2.Split/Zmatch.ExpL.BrainSpanMatch/"
# mkdir -p $OutDir
#parallel -j 20 bash run_Z2.sh {} $Step $ExpMat $MatchDir $OutDir ::: $(seq $Start $Step $End)

# ExpMat="../dat3/HumanCT.Log.SumNorm.Z1.csv"
# MatchDir="/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/MatchGenes.HumanCT/"
# OutDir="/home/jw3514/Work/CellType_Psy/dat2/Z2.Split/SpecMat.Z2.SumNorm.HumanCTMatch/"
# mkdir -p $OutDir
#parallel -j 20 bash run_Z2.sh {} $Step $ExpMat $MatchDir $OutDir ::: $(seq $Start $Step $End)


# ExpMat="../dat2/HumanCT.ExpL.log.csv"
# MatchDir="/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/MatchGenes.HumanCT/"
# OutDir="/home/jw3514/Work/CellType_Psy/dat2/Z2.Split/SpecMat.Z3.HumanCTMatch/"
# mkdir -p $OutDir
# #python $Script --start 0 --step 2 -i $ExpMat -m $MatchDir -o $OutDir
# #parallel -j 20 bash run_Z2.sh {} $Step $ExpMat $MatchDir $OutDir ::: $(seq $Start $Step $End)

# ExpMat="/home/jw3514/Work/CellType_Psy/dat2/Human.CT.AllCell.LogQN.Z1.Entrez.csv"
# MatchDir="/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/MatchGenes.HumanCT/"
# OutDir="/home/jw3514/Work/CellType_Psy/dat2/Z2.Split/SpecMat.logQN.Z2.HumanCTMatch/"
# mkdir -p $OutDir
# #python $Script --start 0 --step 2 -i $ExpMat -m $MatchDir -o $OutDir
# #parallel -j 20 bash run_Z2.sh {} $Step $ExpMat $MatchDir $OutDir ::: $(seq $Start $Step $End)


# ExpMat="/home/jw3514/Work/CellType_Psy/dat2/HumanCT.ExpL.lognorm.csv"
# MatchDir="/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/MatchGenes.HumanCT/"
# OutDir="/home/jw3514/Work/CellType_Psy/dat2/Z2.Split/SpecMat.logQN.Z3.HumanCTMatch/"
# mkdir -p $OutDir
# #python $Script --start 0 --step 2 -i $ExpMat -m $MatchDir -o $OutDir
# #parallel -j 20 bash run_Z2.sh {} $Step $ExpMat $MatchDir $OutDir ::: $(seq $Start $Step $End)

# ExpMat="/home/jw3514/Work/CellType_Psy/dat2/Human.CT.Neuro.LogSubQN.Z1.Entrez.csv"
# MatchDir="/home/jw3514/Work/ASD_Circuits/dat/genes/BrainSpanMatch_uniform_kernal/"
# OutDir="/home/jw3514/Work/CellType_Psy/dat2/Z2.Split/SpecMat.Neuro.logQN.Z2.BrainSpanMatch/"
# #mkdir -p $OutDir
# #nohup parallel -j 19 bash run_Z2.sh {} $Step $ExpMat $MatchDir $OutDir ::: $(seq $Start $Step $End) > test1.log &

# ExpMat="/home/jw3514/Work/CellType_Psy/dat2/Human.CT.Neuro.LogSubQN.Z1.Entrez.csv"
# MatchDir="/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/MatchGenes.HumanCT/"
# OutDir="/home/jw3514/Work/CellType_Psy/dat2/Z2.Split/SpecMat.Neuro.logQN.Z2.HumanCTMatch/"
# #mkdir -p $OutDir
# #nohup parallel -j 19 bash run_Z2.sh {} $Step $ExpMat $MatchDir $OutDir ::: $(seq $Start $Step $End) > test2.log &

# ExpMat="/home/jw3514/Work/CellType_Psy/dat2/Human.CT.NEURO.LogQN.ExpL.Entrez.csv"
# MatchDir="/home/jw3514/Work/ASD_Circuits/dat/genes/BrainSpanMatch_uniform_kernal/"
# OutDir="/home/jw3514/Work/CellType_Psy/dat2/Z2.Split/SpecMat.Neuro.logQN.Z3.BrainSpanMatch/"
# mkdir -p $OutDir
# #nohup parallel -j 19 bash run_Z2.sh {} $Step $ExpMat $MatchDir $OutDir ::: $(seq $Start $Step $End) > test3.log &

# ExpMat="/home/jw3514/Work/CellType_Psy/dat2/Human.CT.NEURO.LogQN.ExpL.Entrez.csv"
# MatchDir="/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/MatchGenes.HumanCT/"
# OutDir="/home/jw3514/Work/CellType_Psy/dat2/Z2.Split/SpecMat.Neuro.logQN.Z3.HumanCTMatch/"
# #mkdir -p $OutDir
# #nohup parallel -j 19 bash run_Z2.sh {} $Step $ExpMat $MatchDir $OutDir ::: $(seq $Start $Step $End) > test4.log &

