#MatchDF="../dat2/ExpMatch/HumanCT.MatchDF.csv"
#OutDir="../dat2/ExpMatch/MatchGenes.HumanCT/"
#mkdir -p $OutDir
#python script_ExpressionMatching.py -i $MatchDF -o $OutDir


#MatchDF="../dat2/ExpMatch/BrainSpan.MatchDF.csv"
#OutDir="../dat2/ExpMatch/MatchGenes.BrainSpan/"
#mkdir -p $OutDir
#python script_ExpressionMatching.py -i $MatchDF -o $OutDir


#MatchDF="/home/jw3514/Work/CellType_Psy/AllenBrainCellAtlas/dat/ABC_TotalExp.Match.csv"
#OutDir="/home/jw3514/Work/CellType_Psy/AllenBrainCellAtlas/dat/MatchGenes/ABC_TotalExp"
#mkdir -p $OutDir
#python script_ExpressionMatching.py -i $MatchDF -o $OutDir


#MatchDF="/home/jw3514/Work/CellType_Psy/AllenBrainCellAtlas/dat/ABC_TotalExp.Match.10xV2.csv"
#OutDir="/home/jw3514/Work/CellType_Psy/AllenBrainCellAtlas/dat/MatchGenes/ABC_TotalExp_10xV2"
#mkdir -p $OutDir
#nohup python script_ExpressionMatching.py -i $MatchDF -o $OutDir &

#MatchDF="/home/jw3514/Work/CellType_Psy/AllenBrainCellAtlas/dat/ABC_TotalExp.Match.10xV3.csv"
#OutDir="/home/jw3514/Work/CellType_Psy/AllenBrainCellAtlas/dat/MatchGenes/ABC_TotalExp_10xV3"
#mkdir -p $OutDir
#nohup python script_ExpressionMatching.py -i $MatchDF -o $OutDir &

MatchDF="/home/jw3514/Work/CellType_Psy/AllenBrainCellAtlas/dat/ABC_TotalExp.Match.CB.csv"
OutDir="/home/jw3514/Work/CellType_Psy/AllenBrainCellAtlas/dat/MatchGenes/ABC_TotalExp_CB/"
mkdir -p $OutDir
nohup python script_ExpressionMatching.py -i $MatchDF -o $OutDir &
