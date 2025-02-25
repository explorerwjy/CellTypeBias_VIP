from ASD_Circuits import *
from scipy.stats import fisher_exact
#from scipy.stats import binom_test
from scipy.stats import hypergeom
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection, multipletests
from scipy import stats

# Load some common used variables
HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()
Anno = pd.read_excel("/home/jw3514/Work/data/HumanBrainCellType/annotation.xlsx", index_col="Cluster")
Anno.drop(Anno.tail(1).index,inplace=True) # drop last n rows
Anno.fillna('', inplace=True)
Anno.index = [int(x) for x in Anno.index.values]

for i, row in Anno.iterrows():
    Class, NT = row["Class auto-annotation"], row["Neurotransmitter auto-annotation"]
    if NT != "":
        Anno.loc[i, "Class auto-annotation"] = "NEUR"
Neur_idx = [int(x) for x in Anno[Anno["Class auto-annotation"]=="NEUR"].index]
NonNeur_idx = [int(x) for x in Anno[Anno["Class auto-annotation"]!="NEUR"].index]

Neurons = sorted(list(set(Anno[Anno["Class auto-annotation"]=="NEUR"]["Supercluster"])))
Not_Neurons = sorted(list(set(Anno[Anno["Class auto-annotation"]!="NEUR"]["Supercluster"])))
Not_Neurons = [x for x in Not_Neurons if x not in Neurons]
ALL_CTs = Neurons + Not_Neurons


def ZscoreConverting_V2(values, mean=np.nan, std=np.nan, low_exp = 0, min_z=-5): 
    """
    Convert values to z-scores with special handling for zeros:
    - Build distribution using only non-zero values
    - Set minimum z-score to min_z (default -5)
    - Set all zero expressions to min_z (default -5)
    
    Args:
        values: Array-like input values
        mean: Optional pre-computed mean
        std: Optional pre-computed standard deviation
        min_z: Minimum z-score (default -5)
    Returns:
        numpy array of z-scores
    """
    # Convert to numpy array and identify non-zero values
    values = np.array(values)
    non_zero_mask = values >= low_exp
    non_zero_values = values[non_zero_mask]
    
    # If no non-zero values, return array of -5
    if len(non_zero_values) == 0:
        return np.full_like(values, min_z)
    
    # Calculate mean and std from non-zero values if not provided
    if mean != mean:  # Check for nan
        mean = np.mean(non_zero_values)
    if std != std:    # Check for nan
        std = np.std(non_zero_values)
        # Handle case where std is 0
        if std == 0:
            std = 1
    
    # Calculate z-scores
    zscores = np.full_like(values, min_z)  # Initialize with -5
    non_zero_zscores = (non_zero_values - mean) / std
    
    # Apply minimum threshold 
    non_zero_zscores = np.maximum(non_zero_zscores, min_z)
    
    # Put non-zero z-scores back in original positions
    zscores[non_zero_mask] = non_zero_zscores
    
    return zscores

def Z1Conversion(ExpMat, outname="test.z1.mat"):
    Z_mat = []
    for g, row in ExpMat.iterrows():
        tmp = ZscoreConverting(row.values)
        Z_mat.append(tmp)
    Z_mat = np.array(Z_mat)
    CT_Z1_DF = pd.DataFrame(data=Z_mat, index=ExpMat.index.values, 
                            columns=ExpMat.columns.values)
    CT_Z1_DF.to_csv(outname)
    return CT_Z1_DF

def Z1Conversion_V2(ExpMat, outname="test.z1.mat", low_exp = 0, min_z=-5):
    """
    Convert expression matrix to z-scores with zero handling
    
    Args:
        ExpMat: pandas DataFrame with genes as rows and cell types as columns
        outname: output file name for saving results
        
    Returns:
        pandas DataFrame with z-scores
    """
    Z_mat = []
    for g, row in ExpMat.iterrows():
        tmp = ZscoreConverting_V2(row.values, min_z=min_z)
        Z_mat.append(tmp)
    
    Z_mat = np.array(Z_mat)
    CT_Z1_DF = pd.DataFrame(data=Z_mat, 
                           index=ExpMat.index.values,
                           columns=ExpMat.columns.values)
    
    CT_Z1_DF.to_csv(outname)
    return CT_Z1_DF


def PlotBiasContrast(X, Y, label1, label2, title=""):
    fig, ax = plt.subplots(dpi=120, figsize=(4.2,4))
    ax.scatter(X, Y, s=15, )

    ax.set_title(title)
    xmin = min(min(X), min(Y), -max(X), -max(Y)) * 1.05
    xmax = max(max(X), max(Y), -min(X), -min(Y)) * 1.05

    ax.plot([xmin,xmax],[xmin,xmax], color="grey", ls="dashed")
    ax.plot([0,0],[xmin,xmax], color="grey", ls="dotted")
    ax.plot([xmin,xmax],[0,0], color="grey", ls="dotted")
    ax.set_xlim((xmin, xmax))
    ax.set_ylim((xmin, xmax))
    r, p = spearmanr(X, Y)
    ax.text(xmin*0.7, xmax*0.7, s="SpearmanR=%.2f\nP= %.1e"%(r,p))
    ax.set_xlabel(label1)
    ax.set_ylabel(label2)
    #print(spearmanr(res["EFFECT_{}".format(name1)].values, res["EFFECT_{}".format(name2)].values))

def PlotBiasContrast_Diff(X, Y, label1, label2, title="", loc=1):
    fig, ax = plt.subplots(dpi=120, figsize=(4.2,4))
    ax.scatter(X, Y, s=15, )

    ax.set_title(title)
    xmin = min(min(X), min(Y), -max(X), -max(Y)) * 1.05
    xmax = max(max(X), max(Y), -min(X), -min(Y)) * 1.05
    
    BiasDiff = np.mean(X) - np.mean(Y)
    t_wil, p_wil = wilcoxon(X, Y)

    ax.plot([xmin,xmax],[xmin,xmax], color="grey", ls="dashed")
    ax.plot([0,0],[xmin,xmax], color="grey", ls="dotted")
    ax.plot([xmin,xmax],[0,0], color="grey", ls="dotted")
    ax.set_xlim((xmin, xmax))
    ax.set_ylim((xmin, xmax))
    r, p = spearmanr(X, Y)
    if loc==1:
        ax.text(xmin*0.7, xmax*0.7, s="BiasDiff=%.2f\nP= %.1e"%(BiasDiff, p_wil))
    elif loc==2:
        ax.text(xmax*0.1, -xmax*0.7, s="BiasDiff=%.2f\nP= %.1e"%(BiasDiff, p_wil))
    ax.set_xlabel(label1)
    ax.set_ylabel(label2)
    
def PlotBiasContrast_Diff(X, Y, label1, label2, title="", loc=1):
    fig, ax = plt.subplots(dpi=150, figsize=(5, 5))
    
    # Improved scatter plot with better color and marker
    ax.scatter(X, Y, s=25, color='dodgerblue', alpha=0.6, edgecolor='k')

    # Title and axis labels with increased font size and weight
    ax.set_title(title, fontsize=16, weight='bold')
    ax.set_xlabel(label1, fontsize=14, weight='bold')
    ax.set_ylabel(label2, fontsize=14, weight='bold')

    # Setting up limits and plot properties
    xmin = min(min(X), min(Y), -max(X), -max(Y)) * 1.05
    xmax = max(max(X), max(Y), -min(X), -min(Y)) * 1.05
    
    BiasDiff = np.mean(X) - np.mean(Y)
    t_wil, p_wil = wilcoxon(X, Y)

    # Plotting the reference lines with explanations
    ax.plot([xmin,xmax],[xmin,xmax], color="grey", ls="dashed", label='Line of Equality')
    ax.plot([0,0],[xmin,xmax], color="grey", ls="dotted", label='X or Y = 0')
    ax.plot([xmin,xmax],[0,0], color="grey", ls="dotted")

    # Setting axis limits
    ax.set_xlim((xmin, xmax))
    ax.set_ylim((xmin, xmax))
    
    # Spearman correlation
    r, p = spearmanr(X, Y)
    
    # Annotate the Bias Difference and Wilcoxon p-value
    if loc == 1:
        ax.text(xmin*0.7, xmax*0.7, s="BiasDiff=%.2f\nP=%.1e"%(BiasDiff, p_wil), fontsize=12, 
                ha='left', va='top', bbox=dict(facecolor='white', alpha=0.5))
    elif loc == 2:
        ax.text(xmax*0.1, -xmax*0.7, s="BiasDiff=%.2f\nP=%.1e"%(BiasDiff, p_wil), fontsize=12, 
                ha='left', va='top', bbox=dict(facecolor='white', alpha=0.5))
    
    # Add a legend to explain the lines
    ax.legend(fontsize=10, loc='upper left')

    # Improve layout and show plot
    plt.tight_layout()
    plt.show()


def CompareSingleCT_New(Bias1, Bias2, CT, name1="1", name2="2", title="", loc=1):
    res = Bias1.join(Bias2, how = 'inner', lsuffix="_{}".format(name1), rsuffix="_{}".format(name2))
    res["Diff"] = res["EFFECT_{}".format(name1)] - res["EFFECT_{}".format(name2)]
    res = res[res["Supercluster_{}".format(name1)]==CT]
    X = res["EFFECT_{}".format(name1)].values
    Y = res["EFFECT_{}".format(name2)].values
    #PlotBiasContrast(X, Y, name1, name2, title=CT)
    PlotBiasContrast_Diff(X, Y, name1, name2, title=CT, loc=loc)

def GetSingeCellBiasCorr(Bias1, Bias2, name1="1", name2="2"):
    res = Bias1.join(Bias2, how = 'inner', lsuffix="_{}".format(name1), rsuffix="_{}".format(name2))
    res["Diff"] = res["EFFECT_{}".format(name1)] - res["EFFECT_{}".format(name2)]
    X = res["EFFECT_{}".format(name1)].values
    Y = res["EFFECT_{}".format(name2)].values
    r, p = spearmanr(X, Y)
    return r,p

def CompareSingleCT2(Bias1, Bias2, CT_idx_Group, names,  name1="1", name2="2"):
    res = Bias1.join(Bias2, how = 'inner', lsuffix="_{}".format(name1), rsuffix="_{}".format(name2))
    res["Diff"] = res["EFFECT_{}".format(name1)] - res["EFFECT_{}".format(name2)]
    xmin = min(min(res["EFFECT_{}".format(name1)].values), min(res["EFFECT_{}".format(name2)].values)) * 1.1
    xmax = max(max(res["EFFECT_{}".format(name1)].values), max(res["EFFECT_{}".format(name2)].values)) * 1.1
    fig, ax = plt.subplots(dpi=200, figsize=(4.2,4))
    for CT_idx, name in zip(CT_idx_Group, names):
        tmp = res.loc[CT_idx,:]
        ax.scatter(tmp["EFFECT_{}".format(name1)].values, tmp["EFFECT_{}".format(name2)].values, s=15, alpha=0.5, label=name)
    #ax.set_title(name)
    ax.plot([xmin,xmax],[xmin,xmax], color="grey", ls="dashed")
    ax.plot([0,0],[xmin,xmax], color="grey", ls="dotted")
    ax.plot([xmin,xmax],[0,0], color="grey", ls="dotted")
    ax.set_xlim((xmin, xmax))
    ax.set_ylim((xmin, xmax))
    r, p = spearmanr(res["EFFECT_{}".format(name1)].values, res["EFFECT_{}".format(name2)].values)
    ax.text(xmin*0.8, xmax*0.8, s="rho=%.2f p= %.1e"%(r,p))
    ax.set_xlabel(name1)
    ax.set_ylabel(name2)
    ax.legend()
    #print(spearmanr(res["EFFECT_{}".format(name1)].values, res["EFFECT_{}".format(name2)].values))

def CompareSingleCT3(DF, Marker, name="",  name1="1", name2="2"):
    tmp = DF
    res = DF
    xmin = min(min(res["EFFECT_{}".format(name1)].values), min(res["EFFECT_{}".format(name2)].values)) * 1.1
    xmax = max(max(res["EFFECT_{}".format(name1)].values), max(res["EFFECT_{}".format(name2)].values)) * 1.1
    fig, ax = plt.subplots(dpi=120, figsize=(4.2,4))
    cmap="Reds"
    ax.scatter(tmp["EFFECT_{}".format(name1)].values, tmp["EFFECT_{}".format(name2)].values, s=30, c=tmp[Marker], cmap=cmap)
    ax.set_title(Marker)
    ax.plot([xmin,xmax],[xmin,xmax], color="grey", ls="dashed")
    ax.plot([0,0],[xmin,xmax], color="grey", ls="dotted")
    ax.plot([xmin,xmax],[0,0], color="grey", ls="dotted")
    ax.set_xlim((xmin, xmax))
    ax.set_ylim((xmin, xmax))
    r, p = spearmanr(res["EFFECT_{}".format(name1)].values, res["EFFECT_{}".format(name2)].values)
    #ax.text(xmin*0.8, xmax*0.8, s="rho=%.2f p= %.1e"%(r,p))
    ax.set_xlabel(name1)
    ax.set_ylabel(name2)
    #print(spearmanr(res["EFFECT_{}".format(name1)].values, res["EFFECT_{}".format(name2)].values))

def CompareCT_x(Bias1, Bias2, name1="1", name2="2", name0 = "",SuperClusters=ALL_CTs, xmin=0, xmax=0, savefig=""):
    res = Bias1.join(Bias2, how = 'inner', lsuffix="_{}".format(name1), rsuffix="_{}".format(name2))
    res["Diff"] = res["EFFECT_{}".format(name1)] - res["EFFECT_{}".format(name2)]
    #res.to_csv("./test/{}_{}_vs_{}.csv".format(name0, name1, name2))
    res = res[res["Supercluster_{}".format(name1)].isin(SuperClusters)]

    print(spearmanr(res["EFFECT_{}".format(name1)].values, res["EFFECT_{}".format(name2)].values))
    #print(pearsonr(res["EFFECT_{}".format(name1)].values, res["EFFECT_{}".format(name2)].values))

    idx = 0
    N_col = 4
    N_rows = math.ceil(len(SuperClusters)/N_col)
    print(len(SuperClusters))
    if len(SuperClusters) > 22:
        fig = plt.figure(figsize=(12, 20), constrained_layout=False)
    else:
        fig = plt.figure(figsize=(12, 16), constrained_layout=False)
    spec = fig.add_gridspec(ncols=N_col, nrows=N_rows)
    if xmin ==0 and xmax==0:
        xmin = min(min(res["EFFECT_{}".format(name1)].values), min(res["EFFECT_{}".format(name2)].values)) * 1.1
        xmax = max(max(res["EFFECT_{}".format(name1)].values), max(res["EFFECT_{}".format(name2)].values)) * 1.1
    for a in range(N_rows):
        for b in range(N_col):
            ax0 = fig.add_subplot(spec[a, b])
            tmp = res[res["Supercluster_{}".format(name1)]==SuperClusters[idx]]
            ax0.scatter(tmp["EFFECT_{}".format(name1)].values, tmp["EFFECT_{}".format(name2)].values, s=10, )
            ax0.set_title(SuperClusters[idx])
            xmin = min(min(tmp["EFFECT_{}".format(name1)].values), min(tmp["EFFECT_{}".format(name2)].values))
            xmax = max(max(tmp["EFFECT_{}".format(name1)].values), max(tmp["EFFECT_{}".format(name2)].values))
            ax0.plot([xmin,xmax],[xmin,xmax], color="grey", ls="dashed")
            #ax0.plot([0,0],[xmin,xmax], color="grey", ls="dotted")
            #ax0.plot([xmin,xmax],[0,0], color="grey", ls="dotted")
            #ax0.set_xlim((xmin, xmax))
            #ax0.set_xlim((xmin, xmax))
            #ax0.set_ylim((0, 5e-3))
            if idx >= len(SuperClusters)-1:
                break
            idx += 1
    plt.tight_layout()
    fig.text(0.5, -0.01, '{} Bias'.format(name1), ha='center', fontsize=20)
    fig.text(-0.02, 0.5, '{} Bias'.format(name2), va='center', rotation='vertical', fontsize=20)
    if savefig != "":
        plt.savefig(savefig, bbox_inches="tight")
    return res, fig

# For each CT dont set x/ylim
def CompareCT_y(Bias1, Bias2, name1="1", name2="2", name0 = "",SuperClusters=ALL_CTs, xmin=0, xmax=0, savefig=""):
    res = Bias1.join(Bias2, how = 'inner', lsuffix="_{}".format(name1), rsuffix="_{}".format(name2))
    res["Diff"] = res["EFFECT_{}".format(name1)] - res["EFFECT_{}".format(name2)]
    #res.to_csv("./test/{}_{}_vs_{}.csv".format(name0, name1, name2))
    res = res[res["Supercluster_{}".format(name1)].isin(SuperClusters)]

    print(spearmanr(res["EFFECT_{}".format(name1)].values, res["EFFECT_{}".format(name2)].values))
    #print(pearsonr(res["EFFECT_{}".format(name1)].values, res["EFFECT_{}".format(name2)].values))

    idx = 0
    N_col = 4
    N_rows = math.ceil(len(SuperClusters)/N_col)
    print(len(SuperClusters))
    if len(SuperClusters) > 22:
        fig = plt.figure(figsize=(12, 20), constrained_layout=False)
    else:
        fig = plt.figure(figsize=(12, 16), constrained_layout=False)
    spec = fig.add_gridspec(ncols=N_col, nrows=N_rows)
    for a in range(N_rows):
        for b in range(N_col):
            ax0 = fig.add_subplot(spec[a, b])
            tmp = res[res["Supercluster_{}".format(name1)]==SuperClusters[idx]]
            ax0.scatter(tmp["EFFECT_{}".format(name1)].values, tmp["EFFECT_{}".format(name2)].values, s=10, )
            ax0.set_title(SuperClusters[idx])
            #ax0.set_ylim((0, 5e-3))
            if idx >= len(SuperClusters)-1:
                break
            idx += 1
    plt.tight_layout()
    fig.text(0.5, -0.01, '{} Bias'.format(name1), ha='center', fontsize=20)
    fig.text(-0.02, 0.5, '{} Bias'.format(name2), va='center', rotation='vertical', fontsize=20)
    if savefig != "":
        plt.savefig(savefig, bbox_inches="tight")
    return res, fig

def CompareSingleCT(Bias1, Bias2, CT, name1="1", name2="2", title=""):
    res = Bias1.join(Bias2, how = 'inner', lsuffix="_{}".format(name1), rsuffix="_{}".format(name2))
    res["Diff"] = res["EFFECT_{}".format(name1)] - res["EFFECT_{}".format(name2)]
    res = res[res["Supercluster_{}".format(name1)]==CT]
    tmp = res
    xmin = min(min(res["EFFECT_{}".format(name1)].values), min(res["EFFECT_{}".format(name2)].values)) * 1.1
    xmax = max(max(res["EFFECT_{}".format(name1)].values), max(res["EFFECT_{}".format(name2)].values)) * 1.1
    fig, ax = plt.subplots(dpi=120, figsize=(4.2,4))
    ax.scatter(tmp["EFFECT_{}".format(name1)].values, tmp["EFFECT_{}".format(name2)].values, s=15, )
    MeanBias1 = tmp["EFFECT_{}".format(name1)].mean()
    MeanBias2 = tmp["EFFECT_{}".format(name2)].mean()

    if title == "":
        ax.set_title(CT)
    else:
        ax.set_title(title)
    ax.plot([xmin,xmax],[xmin,xmax], color="grey", ls="dashed")
    ax.plot([0,0],[xmin,xmax], color="grey", ls="dotted")
    ax.plot([xmin,xmax],[0,0], color="grey", ls="dotted")
    ax.set_xlim((xmin, xmax))
    ax.set_ylim((xmin, xmax))
    r, p = spearmanr(res["EFFECT_{}".format(name1)].values, res["EFFECT_{}".format(name2)].values)
    #ax.text(xmin*0.8, xmax*0.8, s="rho=%.2f p= %.1e \n BiasDiff=%.2f"%(r,p,MeanBias1-MeanBias2))
    ax.set_xlabel(name1)
    ax.set_ylabel(name2)
    #print(spearmanr(res["EFFECT_{}".format(name1)].values, res["EFFECT_{}".format(name2)].values))
    

def CompareSingleCT2(Bias1, Bias2, CT_idx_Group, names,  name1="1", name2="2"):
    res = Bias1.join(Bias2, how = 'inner', lsuffix="_{}".format(name1), rsuffix="_{}".format(name2))
    res["Diff"] = res["EFFECT_{}".format(name1)] - res["EFFECT_{}".format(name2)]
    xmin = min(min(res["EFFECT_{}".format(name1)].values), min(res["EFFECT_{}".format(name2)].values)) * 1.1
    xmax = max(max(res["EFFECT_{}".format(name1)].values), max(res["EFFECT_{}".format(name2)].values)) * 1.1
    fig, ax = plt.subplots(dpi=200, figsize=(4.2,4))
    for CT_idx, name in zip(CT_idx_Group, names):
        tmp = res.loc[CT_idx,:]
        ax.scatter(tmp["EFFECT_{}".format(name1)].values, tmp["EFFECT_{}".format(name2)].values, s=15, alpha=0.5, label=name)
    #ax.set_title(name)
    ax.plot([xmin,xmax],[xmin,xmax], color="grey", ls="dashed")
    ax.plot([0,0],[xmin,xmax], color="grey", ls="dotted")
    ax.plot([xmin,xmax],[0,0], color="grey", ls="dotted")
    ax.set_xlim((xmin, xmax))
    ax.set_ylim((xmin, xmax))
    r, p = spearmanr(res["EFFECT_{}".format(name1)].values, res["EFFECT_{}".format(name2)].values)
    ax.text(xmin*0.8, xmax*0.8, s="rho=%.2f p= %.1e"%(r,p))
    ax.set_xlabel(name1)
    ax.set_ylabel(name2)
    ax.legend()
    #print(spearmanr(res["EFFECT_{}".format(name1)].values, res["EFFECT_{}".format(name2)].values))

def CompareSingleCT3(DF, Marker, name="",  name1="1", name2="2"):
    tmp = DF
    res = DF
    xmin = min(min(res["EFFECT_{}".format(name1)].values), min(res["EFFECT_{}".format(name2)].values)) * 1.1
    xmax = max(max(res["EFFECT_{}".format(name1)].values), max(res["EFFECT_{}".format(name2)].values)) * 1.1
    fig, ax = plt.subplots(dpi=120, figsize=(4.2,4))
    cmap="Reds"
    ax.scatter(tmp["EFFECT_{}".format(name1)].values, tmp["EFFECT_{}".format(name2)].values, s=30, c=tmp[Marker], cmap=cmap)
    ax.set_title(Marker)
    ax.plot([xmin,xmax],[xmin,xmax], color="grey", ls="dashed")
    ax.plot([0,0],[xmin,xmax], color="grey", ls="dotted")
    ax.plot([xmin,xmax],[0,0], color="grey", ls="dotted")
    ax.set_xlim((xmin, xmax))
    ax.set_ylim((xmin, xmax))
    r, p = spearmanr(res["EFFECT_{}".format(name1)].values, res["EFFECT_{}".format(name2)].values)
    #ax.text(xmin*0.8, xmax*0.8, s="rho=%.2f p= %.1e"%(r,p))
    ax.set_xlabel(name1)
    ax.set_ylabel(name2)
    #print(spearmanr(res["EFFECT_{}".format(name1)].values, res["EFFECT_{}".format(name2)].values))

def CompareCT(Bias1, Bias2, name1="1", name2="2", name0 = "",SuperClusters=ALL_CTs, xmin=0, xmax=0, savefig=""):
    res = Bias1.join(Bias2, how = 'inner', lsuffix="_{}".format(name1), rsuffix="_{}".format(name2))
    res["Diff"] = res["EFFECT_{}".format(name1)] - res["EFFECT_{}".format(name2)]
    #res.to_csv("./test/{}_{}_vs_{}.csv".format(name0, name1, name2))
    res = res[res["Supercluster_{}".format(name1)].isin(SuperClusters)]

    print(spearmanr(res["EFFECT_{}".format(name1)].values, res["EFFECT_{}".format(name2)].values))
    #print(pearsonr(res["EFFECT_{}".format(name1)].values, res["EFFECT_{}".format(name2)].values))

    idx = 0
    N_col = 4
    N_rows = math.ceil(len(SuperClusters)/N_col)
    print(len(SuperClusters))
    if len(SuperClusters) > 22:
        fig = plt.figure(figsize=(12, 20), constrained_layout=False, dpi=480)
    else:
        fig = plt.figure(figsize=(12, 16), constrained_layout=False, dpi=480)
    spec = fig.add_gridspec(ncols=N_col, nrows=N_rows)
    if xmin ==0 and xmax==0:
        xmin = min(min(res["EFFECT_{}".format(name1)].values), min(res["EFFECT_{}".format(name2)].values)) * 1.1
        xmax = max(max(res["EFFECT_{}".format(name1)].values), max(res["EFFECT_{}".format(name2)].values)) * 1.1
    for a in range(N_rows):
        for b in range(N_col):
            ax0 = fig.add_subplot(spec[a, b])
            tmp = res[res["Supercluster_{}".format(name1)]==SuperClusters[idx]]
            ax0.scatter(tmp["EFFECT_{}".format(name1)].values, tmp["EFFECT_{}".format(name2)].values, s=10, )
            ax0.set_title(SuperClusters[idx], fontsize=12)
            ax0.plot([xmin,xmax],[xmin,xmax], color="grey", ls="dashed")
            ax0.plot([0,0],[xmin,xmax], color="grey", ls="dotted")
            ax0.plot([xmin,xmax],[0,0], color="grey", ls="dotted")
            ax0.set_xlim((xmin, xmax))
            ax0.set_xlim((xmin, xmax))
            ax0.tick_params(axis='both', which='major', labelsize=12)
            #ax0.set_ylim((0, 5e-3))
            if idx >= len(SuperClusters)-1:
                break
            idx += 1
    fig.text(0.5, -0.01, '{} Bias'.format(name1), ha='center', fontsize=15)
    fig.text(-0.02, 0.5, '{} Bias'.format(name2), va='center', rotation='vertical', fontsize=15)
    fig.suptitle(name0)
    fig.tight_layout()
    if savefig != "":
        plt.savefig(savefig, bbox_inches="tight")
    return res, fig

def AnnotateCTDat(df, Anno):
    for i, row in df.iterrows():
        df.loc[i, "Class"] = Anno.loc[int(i), "Class auto-annotation"]
        df.loc[i, "Supercluster"] = Anno.loc[int(i), "Supercluster"]
        df.loc[i, "Subtype"] = Anno.loc[int(i), "Subtype auto-annotation"]
        df.loc[i, "Neurotransmitter"] = Anno.loc[int(i), "Neurotransmitter auto-annotation"]
        df.loc[i, "Top three regions"] = Anno.loc[int(i), "Top three regions"]
        df.loc[i, "Top three dissections"] = Anno.loc[int(i), "Top three dissections"]
        df.loc[i, "Number of cells"] = Anno.loc[int(i), "Number of cells"] #Top three dissections
        df.loc[i, "Neuropeptide auto-annotation"] = Anno.loc[int(i), "Neuropeptide auto-annotation"] #Top three dissections
    df.index = [int(i) for i in df.index.values]
    return df

def Enrichment(myset, totalset, Ntotal = 20000, method="hypergeom"):
    N = len(myset)
    n = len(myset.intersection(totalset))
    M = Ntotal #20996 #18454 #20000
    m = len(totalset)
    #print(n, N, m, M)
    Rate = (n*M)/(N*m)
    if method == "hypergeom":
        Pvalue = 1 - hypergeom.cdf(n-1, M, m, N)
    elif method == "binom":
        # your code here
        Pvalue = binom_test(n, N, p=m/M)
    elif method == "fisher_exact":
        # your code here
        Odds, Pvalue = fisher_exact([[M-m, m], [N-n, n]], alternative="greater")
    else:
        raise "Unimplemented Error"
    return Rate, Pvalue

def LoadHumanCTAnno():
    Anno = pd.read_excel("/home/jw3514/Work/data/HumanBrainCellType/annotation.xlsx", index_col="Cluster")
    Anno.drop(Anno.tail(1).index,inplace=True) # drop last n rows
    Anno.fillna('', inplace=True)
    Anno.index = [int(x) for x in Anno.index.values]
    return Anno

def AnnoteSubcluster(df, Anno):
    for i, row in df.iterrows():
        idx = int(i.split("-")[0])
        df.loc[i, "Class"] = Anno.loc[idx, "Class"]
        df.loc[i, "Supercluster"] = Anno.loc[idx, "Supercluster"]
        df.loc[i, "Neurotransmitter"] = Anno.loc[idx, "Neurotransmitter"]
        df.loc[i, "Neuropeptide"] = Anno.loc[idx, "Neuropeptide"] #Top three dissections
        df.loc[i, "Top regions"] = Anno.loc[idx, "Top ROIGroupFine"]
        df.loc[i, "Top structures"] = Anno.loc[idx, "Top ROI"]
        df.loc[i, "Number of cells"] = Anno.loc[idx, "Number of cells"] #Top three dissections
    return df



#################################################
# Allen SC data related Functions
#################################################
ABC_ALL_Class = ['01 IT-ET Glut', '02 NP-CT-L6b Glut', '03 OB-CR Glut',
       '04 DG-IMN Glut', '05 OB-IMN GABA', '06 CTX-CGE GABA',
       '07 CTX-MGE GABA', '08 CNU-MGE GABA', '09 CNU-LGE GABA',
       '11 CNU-HYa GABA', '10 LSX GABA', '12 HY GABA', '13 CNU-HYa Glut',
       '14 HY Glut', '15 HY Gnrh1 Glut', '16 HY MM Glut', '17 MH-LH Glut',
       '18 TH Glut', '19 MB Glut', '20 MB GABA', '21 MB Dopa',
       '22 MB-HB Sero', '23 P Glut', '24 MY Glut', '25 Pineal Glut',
       '26 P GABA', '27 MY GABA', '28 CB GABA', '29 CB Glut',
       '30 Astro-Epen', '31 OPC-Oligo', '32 OEC', '33 Vascular',
       '34 Immune']
ABC_nonNEUR = ['30 Astro-Epen', '31 OPC-Oligo', '32 OEC', '33 Vascular', '34 Immune']
def CompareCT_ABC(Bias1, Bias2, name1="1", name2="2", name0 = "",SuperClusters=ABC_ALL_Class, xmin=0, xmax=0, savefig=""):
    res = Bias1.join(Bias2, how = 'inner', lsuffix="_{}".format(name1), rsuffix="_{}".format(name2))
    res["Diff"] = res["EFFECT_{}".format(name1)] - res["EFFECT_{}".format(name2)]
    #res.to_csv("./test/{}_{}_vs_{}.csv".format(name0, name1, name2))
    res = res[res["class_id_label_{}".format(name1)].isin(SuperClusters)]
    
    columns_to_drop_na = ["EFFECT_{}".format(name1), "EFFECT_{}".format(name2)]
    res = res.dropna(subset=columns_to_drop_na)

    print(pearsonr(res["EFFECT_{}".format(name1)].values, res["EFFECT_{}".format(name2)].values))
    print(spearmanr(res["EFFECT_{}".format(name1)].values, res["EFFECT_{}".format(name2)].values))

    idx = 0
    N_col = 4
    N_rows = math.ceil(len(SuperClusters)/N_col)
    #print(len(SuperClusters))
    if len(SuperClusters) > 22:
        fig = plt.figure(figsize=(12, 20), constrained_layout=False)
    else:
        fig = plt.figure(figsize=(12, 16), constrained_layout=False)
    spec = fig.add_gridspec(ncols=N_col, nrows=N_rows)
    if xmin ==0 and xmax==0:
        xmin = min(min(res["EFFECT_{}".format(name1)].values), min(res["EFFECT_{}".format(name2)].values)) * 1.1
        xmax = max(max(res["EFFECT_{}".format(name1)].values), max(res["EFFECT_{}".format(name2)].values)) * 1.1
    for a in range(N_rows):
        for b in range(N_col):
            ax0 = fig.add_subplot(spec[a, b])
            tmp = res[res["class_id_label_{}".format(name1)]==SuperClusters[idx]]
            ax0.scatter(tmp["EFFECT_{}".format(name1)].values, tmp["EFFECT_{}".format(name2)].values, s=10, )
            ax0.set_title(SuperClusters[idx])
            ax0.plot([xmin,xmax],[xmin,xmax], color="grey", ls="dashed")
            ax0.plot([0,0],[xmin,xmax], color="grey", ls="dotted")
            ax0.plot([xmin,xmax],[0,0], color="grey", ls="dotted")
            ax0.set_xlim((xmin, xmax))
            ax0.set_ylim((xmin, xmax))
            if idx >= len(SuperClusters)-1:
                break
            idx += 1
    plt.tight_layout()
    fig.text(0.5, -0.01, '{} Bias'.format(name1), ha='center', fontsize=20)
    fig.text(-0.02, 0.5, '{} Bias'.format(name2), va='center', rotation='vertical', fontsize=20)
    if savefig != "":
        plt.savefig(savefig, bbox_inches="tight")


#################################################
# Go terms related Functions
#################################################
Go2Uniprot = pk.load(open("/home/jw3514/Work/CellType_Psy/dat3/Goterms/Go2Uniprot.pk", 'rb'))
Uniprot2Entrez = pk.load(open("/home/jw3514/Work/CellType_Psy/dat3/Goterms/Uniprot2Entrez.pk", 'rb'))

def GetALLGo(go, GoID):
    Root = go[GoID]
    all_go = Root.get_all_children()
    all_go.add(GoID)
    return all_go
def GetGeneOfGo2(go, GoID, Go2Uniprot=Go2Uniprot):
    goset = GetALLGo(go, GoID)
    Total_Genes = set([])
    for i, tmpgo in enumerate(goset):
        #print(i, tmpgo)
        if tmpgo in Go2Uniprot:
            geneset = set([Uniprot2Entrez.get(x, 0) for x in Go2Uniprot[tmpgo]])
            Total_Genes = Total_Genes.union(geneset)
    return Total_Genes

def CT_Specific_GoTerm_Intersect(CellType, BG_Genes, Z2Bias, CT_Goterm, go, Anno=Anno, topN=100):
    CT_Idx = Anno[Anno["Supercluster"]==CellType].index.values
    tmpmat = Z2Bias.loc[BG_Genes, CT_Idx]
    for g, row in tmpmat.iterrows():
        tmpmat.loc[g, "Mean"] = np.mean(row)
    tmpmat = tmpmat.sort_values("Mean", ascending=False)
    #tmp_genes = tmpmat[tmpmat["Mean"]>1.0].index.values
    #tmp_genes = tmpmat[(tmpmat["Mean"]>=0.5)&(tmpmat['Mean']<=1.0)].index.values
    tmp_genes = tmpmat[(tmpmat["Mean"]>=0.3)&(tmpmat['Mean']<=0.5)].index.values
    print(CellType, len(tmp_genes))
    print([Entrez2Symbol[x] for x in tmp_genes])
    RelatedGos = {}
    CT_Goterm = CT_Goterm.sort_values(CellType, ascending=False)
    for i, row in CT_Goterm.head(topN).iterrows():
        GoID = i
        GoName = row["GoName"]
        Rho = row[CellType]
        GoGenes = GetGeneOfGo2(go, GoID)

        #print(GoGenes)
        InterGenes = GoGenes.intersection(set(tmp_genes))
        InterSymbol = [Entrez2Symbol[x] for x in InterGenes]
        if len(InterGenes) > 0:
            print(GoID, GoName, Rho, InterSymbol)


def SuperClusterBias_BoxPlot(DF, title, NeuroOnly=False, sortby="mean"):
    dat_Z2 = []
    mean_Z2 = []
    if NeuroOnly:
        for _CT in Neurons:
            tmp = DF[DF["Supercluster"] == _CT]
            dat_Z2.append(tmp["EFFECT"].values)
            if sortby == "median":
                mean_Z2.append(np.median(tmp["EFFECT"].values))
            elif sortby == "mean":
                mean_Z2.append(np.mean(tmp["EFFECT"].values))
        mean_Z2 = np.array(mean_Z2)

    else:
        for _CT in ALL_CTs:
            tmp = DF[DF["Supercluster"] == _CT]
            dat_Z2.append(tmp["EFFECT"].values)
            if sortby == "median":
                mean_Z2.append(np.median(tmp["EFFECT"].values))
            elif sortby == "mean":
                mean_Z2.append(np.mean(tmp["EFFECT"].values))
        mean_Z2 = np.array(mean_Z2)


    # Sorting data by the mean values
    sort_idx = np.argsort(mean_Z2)
    show_dat_Z2 = [dat_Z2[x] for x in sort_idx]
    show_CTs = [ALL_CTs[x] for x in sort_idx]

    # Create a figure and axis
    fig, ax = plt.subplots(dpi=120, figsize=(8, 8))

    # Customize the boxplot
    boxprops = dict(linestyle='-', linewidth=1.5, color='blue', facecolor='lightblue')
    medianprops = dict(linestyle='-', linewidth=2.5, color='firebrick')
    meanprops = dict(marker='D', markeredgecolor='black', markerfacecolor='firebrick')

    # Use Seaborn to set the context for publication-quality plots
    sns.set_context("talk", font_scale=1.2)

    # Draw the boxplot
    bp = ax.boxplot(show_dat_Z2, labels=show_CTs, vert=False, patch_artist=True,
                    boxprops=dict(facecolor='lightblue', color='blue'),
                    medianprops=medianprops, meanprops=meanprops, showmeans=True)

    # Add grid lines
    ax.grid(True, linestyle='--', alpha=0.6)

    # Customize the colors of the boxplot elements
    colors = sns.color_palette("Set2", len(show_dat_Z2))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)

    # Add labels and title
    ax.set_xlabel("{} Mutation Bias".format(title), fontsize=14)
    ax.set_ylabel("Supercluster", fontsize=14)
    #ax.set_title("Boxplot of ASD Mutation Bias by Supercluster", fontsize=16)

    # Adjust the tick parameters for better readability
    ax.tick_params(axis='both', which='major', labelsize=12)

    # Show the plot
    plt.tight_layout()
    plt.show()

# MERFISH Related Functions
def FixSubiculum(DF):
    X = DF.loc["Subiculum_dorsal_part"]
    Y = DF.loc["Subiculum_ventral_part"]
    Z = [(X[0]+Y[0])/2, "Hippocampus", 214]
    DF.loc["Subiculum"] = Z
    DF = DF.drop(["Subiculum_dorsal_part", "Subiculum_ventral_part"])
    return DF


def GetExpQ(pvalues):
    sorted_pvalues = np.sort(pvalues)
    expected = np.linspace(0, 1, len(sorted_pvalues), endpoint=False)[1:]
    expected_quantiles = -np.log10(expected)
    observed_quantiles = -np.log10(sorted_pvalues[1:])
    return expected_quantiles, observed_quantiles


def add_class(BiasDF, ClusterAnn):
    for cluster, row in BiasDF.iterrows():
        BiasDF.loc[cluster, "class_id_label"] = ClusterAnn.loc[cluster, "class_id_label"]
        BiasDF.loc[cluster, "CCF_broad.freq"] = ClusterAnn.loc[cluster, "CCF_broad.freq"]
        BiasDF.loc[cluster, "CCF_acronym.freq"] = ClusterAnn.loc[cluster, "CCF_acronym.freq"]
        BiasDF.loc[cluster, "v3.size"] = ClusterAnn.loc[cluster, "v3.size"]
        BiasDF.loc[cluster, "v2.size"] = ClusterAnn.loc[cluster, "v2.size"]
    return BiasDF


####################################################
# Specificity Score Validation Functions
####################################################
def MergeBiasDF(Bias1, Bias2, name1="1", name2="2"):
    res = Bias1.join(Bias2, how = 'inner', lsuffix="_{}".format(name1), rsuffix="_{}".format(name2))
    res["Diff"] = res["EFFECT_{}".format(name1)] - res["EFFECT_{}".format(name2)]
    return res 
    #res.to_csv("./test/{}_{}_vs_{}.csv".format(name0, name1, name2))

def PlotBiasContrast_v2(MergeDF, name1, name2, dataset="Human", title=""):
    if dataset == "HumanCT":
        NEUR = MergeDF[MergeDF["Supercluster_{}".format(name1)].isin(Neurons)]
        NonNEUR = MergeDF[~MergeDF["Supercluster_{}".format(name1)].isin(Neurons)]
        X_NEUR, Y_NEUR = NEUR["EFFECT_{}".format(name1)], NEUR["EFFECT_{}".format(name2)]
        X_NonNEUR, Y_NonNEUR = NonNEUR["EFFECT_{}".format(name1)], NonNEUR["EFFECT_{}".format(name2)]
        fig, ax = plt.subplots(dpi=300, figsize=(5, 5))

        ax.scatter(X_NEUR, Y_NEUR, s=40, color="blue", edgecolor='black', alpha=0.7)
        ax.scatter(X_NonNEUR, Y_NonNEUR, s=40, color="red", edgecolor='black', alpha=0.7)
        
        ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
        
        # Spearman correlation
        r_Neur, p_Neur = spearmanr(X_NEUR, Y_NEUR)
        r_nonNeur, p_nonNeur = spearmanr(X_NonNEUR, Y_NonNEUR)
        
        xmin = np.min(X_NEUR)
        xmax = np.max(X_NEUR)
        ymax = np.max(Y_NEUR)
        # Adjust text position to avoid overlap
        text_x = xmin * 0.9 if r_Neur < 0 else xmin * 0.7
        text_y = ymax * 0.7
        
        ax.text(text_x, text_y, s=f"R_NEUR = {r_Neur:.2f}\nR_NonNEUR = {r_nonNeur:.2f}",
                fontsize=12, bbox=dict(facecolor='white', alpha=0.5, edgecolor='black'))
        
        ax.set_xlabel(name1, fontsize=12, fontweight='bold')
        ax.set_ylabel(name2, fontsize=12, fontweight='bold')
        
        # Style adjustments
        ax.grid(True, linestyle='--', alpha=0.6)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.xticks(fontsize=10, fontweight='bold')
        plt.yticks(fontsize=10, fontweight='bold')
        
        plt.tight_layout()
        plt.show()
        return
    elif dataset == "MouseCT":
        NEUR = MergeDF[~MergeDF["class_id_label_{}".format(name1)].isin(ABC_nonNEUR)]
        NonNEUR = MergeDF[MergeDF["class_id_label_{}".format(name1)].isin(ABC_nonNEUR)]
        X_NEUR, Y_NEUR = NEUR["EFFECT_{}".format(name1)], NEUR["EFFECT_{}".format(name2)]
        X_NonNEUR, Y_NonNEUR = NonNEUR["EFFECT_{}".format(name1)], NonNEUR["EFFECT_{}".format(name2)]
        fig, ax = plt.subplots(dpi=300, figsize=(5, 5))

        ax.scatter(X_NEUR, Y_NEUR, s=40, color="blue", edgecolor='black', alpha=0.7)
        ax.scatter(X_NonNEUR, Y_NonNEUR, s=40, color="red", edgecolor='black', alpha=0.7)
        
        ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
        
        # Spearman correlation
        r_Neur, p_Neur = spearmanr(X_NEUR, Y_NEUR)
        r_nonNeur, p_nonNeur = spearmanr(X_NonNEUR, Y_NonNEUR)
        
        xmin = np.min(X_NEUR)
        xmax = np.max(X_NEUR)
        ymax = np.max(Y_NEUR)
        # Adjust text position to avoid overlap
        text_x = xmin * 0.9 if r_Neur < 0 else xmin * 0.7
        text_y = ymax * 0.7
        
        ax.text(text_x, text_y, s=f"R_NEUR = {r_Neur:.2f}\nR_NonNEUR = {r_nonNeur:.2f}",
                fontsize=12, bbox=dict(facecolor='white', alpha=0.5, edgecolor='black'))
        
        ax.set_xlabel(name1, fontsize=12, fontweight='bold')
        ax.set_ylabel(name2, fontsize=12, fontweight='bold')
        
        # Style adjustments
        ax.grid(True, linestyle='--', alpha=0.6)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.xticks(fontsize=10, fontweight='bold')
        plt.yticks(fontsize=10, fontweight='bold')
        
        plt.tight_layout()
        plt.show()
        return

# Phenotype vs IQ Functions

def linear_fit(biases, IQs, alpha=0.05, fixed_intercept=None):
    if fixed_intercept is None:
        model = sm.OLS(IQs, sm.add_constant(biases))
        results = model.fit()
        intercept = results.params[0]
    else:
        # Fit model without intercept and add fixed intercept
        model = sm.OLS(IQs - fixed_intercept, biases)
        results = model.fit()
        intercept = fixed_intercept
    
    beta = results.params[0] if fixed_intercept is not None else results.params[1]
    # get CI of beta
    ci = results.conf_int(alpha=alpha)
    ci_low = ci[0][0] if fixed_intercept is not None else ci[1][0]
    ci_high = ci[0][1] if fixed_intercept is not None else ci[1][1]
    r_value = results.rsquared
    p_value = results.pvalues[0] if fixed_intercept is not None else results.pvalues[1]
    std_err = results.bse[0] if fixed_intercept is not None else results.bse[1]
    pho, p = spearmanr(biases, IQs)
    
    return intercept, beta, ci_low, ci_high, r_value, p_value, std_err, pho, p


def Plot_Bias_vs_IQ(STR, Mut_n_IQ_conf, BiasMat):
    biases, IQs = BiasVsPheno(Mut_n_IQ_conf, BiasMat , STR, 'XX')
    pho, p = spearmanr(biases, IQs)
    # Create the scatter plot
    plt.figure(dpi=150, figsize=(5, 4))

    # Plot data points
    plt.scatter(biases, IQs, s=50, color="#2c7bb6", edgecolor="black", alpha=0.8, zorder=10)

    # Fit and plot the trend line
    b, a = np.polyfit(biases, IQs, deg=1)
    xseq = np.linspace(min(biases), max(biases), num=100)
    plt.plot(xseq, a + b * xseq, color="#d7191c", lw=2.5, linestyle='--', zorder=5)
    _SuperCluster = Anno.loc[STR, "Supercluster"]
    # Add title with improved formatting
    plt.title(f'{_SuperCluster} - {STR} \nSpearman ρ = {pho:.2f}, p = {p:.2e}', fontsize=14, fontweight='bold')

    # Labeling axes
    plt.xlabel("Cell Type Bias", fontsize=12, fontweight='bold')

    plt.ylabel("Full Scale IQ", fontsize=12, fontweight='bold')

    # Grid lines for better readability
    plt.grid(True, linestyle='--', alpha=0.5)

    # Adjust tick parameters
    plt.xticks(fontsize=10, fontweight='bold')
    plt.yticks(fontsize=10, fontweight='bold')

    # Tight layout for better spacing
    plt.tight_layout()

    # Show the plot
    plt.show()

def BiasVsPheno(MutPhenoDF, BiasMat, STR, label):
    biases = []
    IQs = []
    for i, row in MutPhenoDF.iterrows():
        if label == 'label':
            gene = row["HGNC"]
        else:
            gene = int(row["Entrez"])
        if gene in BiasMat.index.values:
            bias = BiasMat.loc[gene, STR]
            if bias == bias:
                IQ = row["IQ"]
                biases.append(bias)
                IQs.append(IQ)
    return np.array(biases), np.array(IQs)

def ADJ_P(DF, p_col="p_value"):
    #_, q = fdrcorrection(DF["p"].values)
    #_, q, alphacSidak, alphacBonf= multipletests(DF["p"].values, method='fdr_by')
    DF = DF.sort_values(by=p_col, ascending=True)
    _, q, alphacSidak, alphacBonf= multipletests(DF[p_col].values, method='fdr_bh')
    DF[f'{p_col}_adj'] = q
    return DF

def Plot_Bias_vs_IQ_MoustCT(STR, Mut_n_IQ_conf, HCT_Z2_MAT_HCT, ax=None):
    biases, IQs = BiasVsPheno(Mut_n_IQ_conf, HCT_Z2_MAT_HCT , STR, 'XX')
    pho, p = spearmanr(biases, IQs)
    
    if ax is None:
        fig, ax = plt.subplots(dpi=150, figsize=(5, 4))

    # Plot data points
    ax.scatter(biases, IQs, s=50, color="#2c7bb6", edgecolor="black", alpha=0.8, zorder=10)

    # Fit and plot the trend line
    b, a = np.polyfit(biases, IQs, deg=1)
    xseq = np.linspace(min(biases), max(biases), num=100)
    ax.plot(xseq, a + b * xseq, color="#d7191c", lw=2.5, linestyle='--', zorder=5)

  

    # Compute average bias for IQ > 70 and IQ < 70
    high_iq_bias = np.mean([bias for bias, iq in zip(biases, IQs) if iq > 70])
    low_iq_bias = np.mean([bias for bias, iq in zip(biases, IQs) if iq <= 70])
    bias_diff = high_iq_bias - low_iq_bias

    # Add title with improved formatting and average bias information
    # ax.set_title(f'{STR}\nSpearman ρ = {pho:.2f}, p = {p:.2e}\n'
    #              f'Avg Bias (IQ>70): {high_iq_bias:.2f}, (IQ≤70): {low_iq_bias:.2f}, Diff: {bias_diff:.2f}', 
    #              fontsize=12, fontweight='bold')

    ax.set_title(f'{STR}\nSpearman ρ = {pho:.2f}, p = {p:.2e} PBS = {b:.2f}', 
                 fontsize=12, fontweight='bold')
  # Add line at IQ 70
    #ax.axhline(y=70, color='green', linestyle=':', linewidth=2)
    # plot HIQ, LIQ and bias diff as x-axis on y=70
    #ax.plot([high_iq_bias, low_iq_bias], [70, 70], color='black', linestyle='-', linewidth=5)
    # Labeling axes
    ax.set_xlabel("Cell Type Bias", fontsize=12, fontweight='bold')
    ax.set_ylabel("Full Scale IQ", fontsize=12, fontweight='bold')

    # Grid lines for better readability
    ax.grid(True, linestyle='--', alpha=0.5)

    # Adjust tick parameters
    ax.tick_params(axis='both', which='major', labelsize=10)

    # Tight layout for better spacing
    plt.tight_layout()

    return ax

def Plot_Bias_vs_IQ_MoustCT_forPaper(STR, Mut_n_IQ_conf, HCT_Z2_MAT_HCT, fixed_intercept=None, ax=None):
    biases, IQs = BiasVsPheno(Mut_n_IQ_conf, HCT_Z2_MAT_HCT , STR, 'XX')
    intercept, beta, ci_low, ci_high, r_value, p_value, std_err, pho, p = linear_fit(biases, IQs, alpha=0.05, fixed_intercept=fixed_intercept)
    #pho, p = spearmanr(biases, IQs)
    
    if ax is None:
        fig, ax = plt.subplots(dpi=150, figsize=(5, 4))

    # Plot data points
    ax.scatter(biases, IQs, s=50, color="#2c7bb6", edgecolor="black", alpha=0.8, zorder=10)

    # Fit and plot the trend line
    #b, a = np.polyfit(biases, IQs, deg=1)
    xseq = np.linspace(min(biases), max(biases), num=100)
    ax.plot(xseq, intercept + beta * xseq, color="#d7191c", lw=2.5, linestyle='--', zorder=5)

    # Compute average bias for IQ > 70 and IQ < 70
    high_iq_bias = np.mean([bias for bias, iq in zip(biases, IQs) if iq > 70])
    low_iq_bias = np.mean([bias for bias, iq in zip(biases, IQs) if iq <= 70])
    bias_diff = high_iq_bias - low_iq_bias


    #ax.set_title(f'{STR}\nSpearman ρ = {pho:.2f}, p = {p:.2e} \nPBS = {beta:.2f} PBS_p_value = {p_value:.2e}', 
    #            fontsize=12, fontweight='bold')
    
    #ax.set_title(f'{STR}\nPBS = {beta:.2f} p_value = {p_value:.1e}', 
    #             fontsize=12, fontweight='bold')
    ax.set_title(f'{STR}\nPBS = {beta:.2f} p_value = {p_value:.0e}', 
                 fontsize=12, fontweight='bold')
  # Add line at IQ 70
    #ax.axhline(y=70, color='green', linestyle=':', linewidth=2)
    # plot HIQ, LIQ and bias diff as x-axis on y=70
    #ax.plot([high_iq_bias, low_iq_bias], [70, 70], color='black', linestyle='-', linewidth=5)
    # Labeling axes
    ax.set_xlabel("Cell Type Bias", fontsize=12, fontweight='bold')
    ax.set_ylabel("Full Scale IQ", fontsize=12, fontweight='bold')

    # Grid lines for better readability
    ax.grid(True, linestyle='--', alpha=0.5)

    # Adjust tick parameters
    ax.tick_params(axis='both', which='major', labelsize=10)

    # Tight layout for better spacing
    plt.tight_layout()

    return ax

def Plot_Bias_vs_IQ_HumanCT(STR, Mut_n_IQ_conf, HCT_Z2_MAT_HCT, ax=None, Pval=None):
    biases, IQs = BiasVsPheno(Mut_n_IQ_conf, HCT_Z2_MAT_HCT , STR, 'XX')
    pho, p = spearmanr(biases, IQs)
    if Pval is None:
        Pval = p
    if ax is None:
        fig, ax = plt.subplots(dpi=150, figsize=(5, 4))
    ax.scatter(biases, IQs, s=50, color="#2c7bb6", edgecolor="black", alpha=0.8, zorder=10)

    b, a = np.polyfit(biases, IQs, deg=1)
    xseq = np.linspace(min(biases), max(biases), num=100)
    ax.plot(xseq, a + b * xseq, color="#d7191c", lw=2.5, linestyle='--', zorder=5)
    # Compute average bias for IQ > 70 and IQ < 70
    high_iq_bias = np.mean([bias for bias, iq in zip(biases, IQs) if iq > 70])
    low_iq_bias = np.mean([bias for bias, iq in zip(biases, IQs) if iq <= 70])
    bias_diff = high_iq_bias - low_iq_bias

    _SuperCluster = Anno.loc[STR, "Supercluster"]

    ax.set_title(f'{_SuperCluster} - {STR}', 
            fontsize=25, fontweight='normal')
    ax.text(0.60, 0.85, s=f'PBS = {b:.2f}\np = {Pval:.1e}',
            fontsize=22.5, ha='left', va='top', transform=ax.transAxes)

    ax.set_xlabel("Cell Type Bias", fontsize=25, fontweight='normal')
    ax.set_ylabel("Full Scale IQ", fontsize=25, fontweight='normal')

    # Grid lines for better readability
    ax.grid(True, linestyle='--', alpha=0.5)

    # Adjust tick parameters
    ax.tick_params(axis='both', which='major', labelsize=15)

    return ax



# Whole Genome Regression 
def plot_cluster_correlation(cluster, SCZMutDF, specificity_scores, eff_label = "LGD_OR", plot=False):
    """
    Plot correlation between cluster bias and LGD odds ratio.
    
    Args:
        cluster (int): Cluster number to analyze
        SCZMutDF (pd.DataFrame): DataFrame containing mutation data
        specificity_scores (pd.DataFrame): Matrix of bias/specificity scores
    """
    entrez_list = SCZMutDF.index.values
    Zscore_list = SCZMutDF[eff_label].values
    Bias_list = specificity_scores.loc[entrez_list, cluster]
    valid_mask = ~np.isnan(Zscore_list) & ~np.isnan(Bias_list)

    # Calculate correlations
    spearman_corr, spearman_p = stats.spearmanr(Zscore_list[valid_mask], Bias_list[valid_mask])
    pearson_corr, pearson_p = stats.pearsonr(Zscore_list[valid_mask], Bias_list[valid_mask])
    #print(f"Spearman correlation: {spearman_corr}")
    #print(f"Pearson correlation: {pearson_corr}")

    # Clean data
    Zscore_list_clean = Zscore_list[valid_mask]
    Bias_list_clean = Bias_list[valid_mask]

    # Fit linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(Bias_list_clean, Zscore_list_clean, alternative="greater")
    
    if plot:
        # Create scatter plot
        plt.scatter(Bias_list_clean, Zscore_list_clean, s=1)

        # Add regression line
        x_range = np.array([min(Bias_list_clean), max(Bias_list_clean)])
        plt.plot(x_range, slope * x_range + intercept, 'r', 
                label=f'β = {slope:.2f}\nR² = {r_value**2:.3f}\np = {p_value:.2e}')

        plt.xlabel("Bias")
        plt.ylabel("Z-score") 
        plt.legend()
        plt.show()
    return spearman_corr, spearman_p, pearson_corr, pearson_p, slope,std_err, r_value, p_value, 
# Example usage:

def cell_type_bias_Linear_fit(specificity_scores, Spark_Meta_test, Anno, eff_label="ZSTAT"):

    intercept, beta, ci_low, ci_high, r_value, p_value, std_err, pho, p = linear_fit(specificity_scores, Spark_Meta_test, Anno, eff_label=eff_label)
    return beta, r_value, std_err, p_value

def calculate_cluster_correlations(specificity_scores, Spark_Meta_test, Anno, eff_label="AutismMerged_LoF"):
    # Create lists to store results
    clusters = []
    spearman_correlations = []
    spearman_pvalues = []
    pearson_correlations = []
    pearson_pvalues = []
    superclusters = []
    slope_values = []
    std_err_values = []
    r_value_values = []
    p_value_values = []


    for cluster in specificity_scores.columns.values:
        #spearman_corr, spearman_p, pearson_corr, pearson_p = plot_cluster_correlation(cluster, Spark_Meta_test, specificity_scores, eff_label=eff_label)
        spearman_corr, spearman_p, pearson_corr, pearson_p, slope, std_err, r_value, p_value = plot_cluster_correlation(cluster, Spark_Meta_test, specificity_scores, eff_label=eff_label)
        # Store results
        clusters.append(cluster)
        spearman_correlations.append(spearman_corr)
        spearman_pvalues.append(spearman_p)
        #pearson_correlations.append(pearson_corr)
        #pearson_pvalues.append(pearson_p)
        slope_values.append(slope)
        std_err_values.append(std_err)
        r_value_values.append(r_value)
        p_value_values.append(p_value)
        superclusters.append(Anno.loc[cluster, "Supercluster"])

    # Create DataFrame with results
    corr_df_ASD = pd.DataFrame({
        'Cluster': clusters,
        'SuperCluster': superclusters,
        'Spearman_Correlation': spearman_correlations,
        'Spearman_P_value': spearman_pvalues,
        'Slope': slope_values,
        'Std_err': std_err_values,
        'R_value': r_value_values,
        'P_value': p_value_values
    })
    corr_df_ASD = corr_df_ASD.sort_values(by="Spearman_P_value", ascending=True)
    return corr_df_ASD

def plot_supercluster_correlations(corr_df, title=""):
    """
    Create a box plot showing distribution of Spearman correlations by SuperCluster.
    
    Args:
        corr_df: DataFrame containing 'SuperCluster' and 'Spearman_Correlation' columns
    """
    # Calculate mean correlation per SuperCluster for sorting
    mean_corr = corr_df.groupby('SuperCluster')['Spearman_Correlation'].mean().sort_values(ascending=False)
    order = mean_corr.index

    # Create a box plot grouped by SuperCluster
    plt.figure(figsize=(12, 10))
    sns.boxplot(data=corr_df, y='SuperCluster', x='Spearman_Correlation', order=order)
    plt.title(title)
    plt.tight_layout()
    plt.show()