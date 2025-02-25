#!/bin/bash

#run simulation in parallel of 100X100 each

GeneProbFil="/home/jw3514/Work/ASD_Circuits/dat/genes/May14_Gene_n_Prob.LGD.Dmis.csv"
Z2Mat="../dat/HumanCellType.AllCell.HCT.Z2bias.entrez.csv"

#WeightFil="/home/jw3514/Work/ASD_Circuits/dat/Unionize_bias/SCZ.top200.gw.csv"
#WeightFil="/home/jw3514/Work/CellType_Psy/dat3/SCZ.top200.gw.csv"
#Location="../dat2/Bias/SubSampled.Sibling.SCZ200.Weights/"
#mkdir -p $Location
#parallel -j 10 bash run_cont_rand_weight_transfer.sh -a {} -w $WeightFil -l $Location -p $GeneProbFil -i $Z2Mat ::: $(seq 1 10)

#WeightFil="/home/jw3514/Work/ASD_Circuits/dat/Unionize_bias/Spark_Meta_EWS.GeneWeight.csv"
WeightFil="/home/jw3514/Work/CellType_Psy/dat3/ASD.top200.gw.csv"
Location="../dat2/Bias/SubSampled.Sibling.ASD200.Weights/"
mkdir -p $Location
parallel -j 10 bash run_cont_rand_weight_transfer.sh -a {} -w $WeightFil -l $Location -p $GeneProbFil -i $Z2Mat ::: $(seq 1 10)


#WeightFil="/home/jw3514/Work/CellType_Psy/dat/Bipolar.top200.gw.csv"
#Location="../dat2/Bias/SubSampled.Sibling.BP200.Weights/"
#mkdir -p $Location
#parallel -j 10 bash run_cont_rand_weight_transfer.sh -a {} -w $WeightFil -l $Location -p $GeneProbFil -i $Z2Mat ::: $(seq 1 10)


#WeightFil="/home/jw3514/Work/CellType_Psy/dat3/DDD.hc.asd.scz.exclude.gw.csv"
#Location="../dat2/Bias/SubSampled.Sibling.NDD.Weights/"
#mkdir -p $Location
#parallel -j 10 bash run_cont_rand_weight_transfer.sh -a {} -w $WeightFil -l $Location -p $GeneProbFil -i $Z2Mat ::: $(seq 1 10)


