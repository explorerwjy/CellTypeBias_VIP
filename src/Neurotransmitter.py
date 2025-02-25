# Neurotransmitter.py
#Neurotransmitter Systems Dictionary
import sys
sys.path.insert(1, '/home/jw3514/Work/CellType_Psy/src')
from CellType_PSY import *

# Dictionary of neurotransmitter systems with associated genes
neurotransmitter_systems = {
    "glutamate": {
        "type": "excitatory",
        "synthesis_transport": {
            "vesicular_transporters": ["SLC17A6", "SLC17A7"],  # VGLUT2, VGLUT1
            "synthesis": ["GLS", "GLS2"],
            "transporters": ["SLC1A1", "SLC1A2", "SLC1A3", "SLC1A4", "SLC1A5", "SLC1A6", "SLC1A7"]
        },
        "receptors": {
            "ionotropic": {
                "NMDA": ["GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D"],
                "AMPA": ["GRIA1", "GRIA2", "GRIA3", "GRIA4"],
                "kainate": ["GRIK1", "GRIK2", "GRIK3", "GRIK4", "GRIK5"]
            },
            "metabotropic": {
                "mGluR": ["GRM1", "GRM2", "GRM3", "GRM4", "GRM5", "GRM6", "GRM7", "GRM8"]
            }
        }
    },
    
    "GABA": {
        "type": "inhibitory",
        "synthesis_transport": {
            "synthesis": ["GAD1", "GAD2"],  # GAD67, GAD65
            "vesicular_transporter": ["SLC32A1"],  # VGAT
            "transporter": ["SLC6A1"]  # GAT1
        },
        "receptors": {
            "GABA_A": {
                "alpha": ["GABRA1", "GABRA2", "GABRA3", "GABRA4", "GABRA5", "GABRA6"],
                "beta": ["GABRB1", "GABRB2", "GABRB3"],
                "gamma": ["GABRG1", "GABRG2", "GABRG3"],
                "delta": ["GABRD"]
            },
            "GABA_B": ["GABBR1", "GABBR2"]
        }
    },
    
    "dopamine": {
        "type": "modulatory",
        "synthesis_transport": {
            "synthesis": ["TH", "DDC"],
            "transporter": ["SLC6A3"],  # DAT
            "vesicular_transporter": ["SLC18A2"]  # VMAT2
        },
        "receptors": {
            "D1_like": ["DRD1", "DRD5"],
            "D2_like": ["DRD2", "DRD3", "DRD4"]
        }
    },
    
    "serotonin": {
        "type": "modulatory",
        "synthesis_transport": {
            "synthesis": ["TPH1", "TPH2", "DDC"],
            "transporter": ["SLC6A4"],  # SERT
            "vesicular_transporter": ["SLC18A2"]  # VMAT2
        },
        "receptors": {
            "5HT1": ["HTR1A", "HTR1B", "HTR1D", "HTR1E", "HTR1F"],
            "5HT2": ["HTR2A", "HTR2B", "HTR2C"],
            "5HT3": ["HTR3A", "HTR3B", "HTR3C", "HTR3D", "HTR3E"],
            "others": ["HTR4", "HTR5A", "HTR6", "HTR7"]
        }
    },
    
    "acetylcholine": {
        "type": "modulatory",
        "synthesis_transport": {
            "synthesis": ["CHAT"],
            "transporter": ["SLC5A7"],  # CHT1
            "vesicular_transporter": ["SLC18A3"],  # VAChT
            "degradation": ["ACHE"]
        },
        "receptors": {
            "nicotinic": {
                "alpha": ["CHRNA1", "CHRNA2", "CHRNA3", "CHRNA4", "CHRNA5", 
                         "CHRNA6", "CHRNA7", "CHRNA8", "CHRNA9", "CHRNA10"],
                "beta": ["CHRNB1", "CHRNB2", "CHRNB3", "CHRNB4"],
                "other": ["CHRND", "CHRNE", "CHRNG"]
            },
            "muscarinic": ["CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5"]
        }
    },
    
    "norepinephrine": {
        "type": "modulatory",
        "synthesis_transport": {
            "synthesis": ["TH", "DDC", "DBH"],
            "transporter": ["SLC6A2"],  # NET
            "vesicular_transporter": ["SLC18A2"]  # VMAT2
        },
        "receptors": {
            "alpha1": ["ADRA1A", "ADRA1B", "ADRA1D"],
            "alpha2": ["ADRA2A", "ADRA2B", "ADRA2C"],
            "beta": ["ADRB1", "ADRB2", "ADRB3"]
        }
    },
    
    "regulatory_genes": {
        "degradation": ["MAOA", "MAOB", "COMT"],
        "synthesis_regulation": ["DBH", "PNMT"]
    }
}
neurotransmitter_systems["oxytocin"] = {
    "type": "neuromodulator",
    "synthesis_transport": {
        "synthesis": ["OXT"],
        "secretion_regulation": ["CD38", "CD157"],
        "degradation": ["LNPEP"]
    },
    "receptors": {
        "primary_receptor": ["OXTR"]
    }
}

gene_to_entrez = {
    # Glutamate related
    "SLC17A6": 57084, "SLC17A7": 57030,  # VGLUTs
    "GLS": 2744, "GLS2": 27165,
    "SLC1A1": 6505, "SLC1A2": 6506, "SLC1A3": 6507, 
    "SLC1A4": 6509, "SLC1A5": 6510, "SLC1A6": 6511, "SLC1A7": 6512,
    "GRIN1": 2902, "GRIN2A": 2903, "GRIN2B": 2904, "GRIN2C": 2905, "GRIN2D": 2906,
    "GRIA1": 2890, "GRIA2": 2891, "GRIA3": 2892, "GRIA4": 2893,
    "GRIK1": 2897, "GRIK2": 2898, "GRIK3": 2899, "GRIK4": 2900, "GRIK5": 2901,
    "GRM1": 2911, "GRM2": 2912, "GRM3": 2913, "GRM4": 2914,
    "GRM5": 2915, "GRM6": 2916, "GRM7": 2917, "GRM8": 2918,
    
    # GABA related
    "GAD1": 2571, "GAD2": 2572,  # GAD67, GAD65
    "SLC32A1": 140679,  # VGAT
    "SLC6A1": 6529,  # GAT1
    "GABRA1": 2554, "GABRA2": 2555, "GABRA3": 2556,
    "GABRA4": 2557, "GABRA5": 2558, "GABRA6": 2559,
    "GABRB1": 2560, "GABRB2": 2561, "GABRB3": 2562,
    "GABRG1": 2565, "GABRG2": 2566, "GABRG3": 2567,
    "GABRD": 2563,
    "GABBR1": 2550, "GABBR2": 9568,
    
    # Dopamine related
    "TH": 7054,
    "DDC": 1644,
    "SLC6A3": 6531,  # DAT
    "SLC18A2": 6571,  # VMAT2
    "DRD1": 1812, "DRD2": 1813, "DRD3": 1814,
    "DRD4": 1815, "DRD5": 1816,
    
    # Serotonin related
    "TPH1": 7166, "TPH2": 121278,
    "SLC6A4": 6532,  # SERT
    "HTR1A": 3350, "HTR1B": 3351, "HTR1D": 3352,
    "HTR1E": 3354, "HTR1F": 3355,
    "HTR2A": 3356, "HTR2B": 3357, "HTR2C": 3358,
    "HTR3A": 3359, "HTR3B": 9177, "HTR3C": 170572,
    "HTR3D": 200909, "HTR3E": 285242,
    "HTR4": 3360, "HTR5A": 3361, "HTR6": 3362, "HTR7": 3363,
    
    # Acetylcholine related
    "CHAT": 1103,
    "SLC5A7": 60482,  # CHT1
    "SLC18A3": 6572,  # VAChT
    "ACHE": 43,
    "CHRNA1": 1134, "CHRNA2": 1135, "CHRNA3": 1136,
    "CHRNA4": 1137, "CHRNA5": 1138, "CHRNA6": 8973,
    "CHRNA7": 1139, "CHRNA8": 57053, "CHRNA9": 55584,
    "CHRNA10": 57053,
    "CHRNB1": 1140, "CHRNB2": 1141, "CHRNB3": 1142, "CHRNB4": 1143,
    "CHRND": 1144, "CHRNE": 1145, "CHRNG": 1146,
    "CHRM1": 1128, "CHRM2": 1129, "CHRM3": 1131,
    "CHRM4": 1132, "CHRM5": 1133,
    
    # Norepinephrine related
    "DBH": 1621,
    "SLC6A2": 6530,  # NET
    "ADRA1A": 148, "ADRA1B": 147, "ADRA1D": 146,
    "ADRA2A": 150, "ADRA2B": 151, "ADRA2C": 152,
    "ADRB1": 153, "ADRB2": 154, "ADRB3": 155,
    
    # Regulatory genes
    "MAOA": 4128, "MAOB": 4129,
    "COMT": 1312,
    "PNMT": 5409
}
gene_to_entrez.update({
    # Oxytocin synthesis and processing
    "OXT": 5020,      # Oxytocin/neurophysin I prepropeptide
    "OXTR": 5021,     # Oxytocin receptor
    "OXTL": 5022,     # Oxytocin-like gene
    "CD38": 952,      # CD38 molecule - involved in oxytocin secretion
    "LNPEP": 4012,    # Leucyl/cystinyl aminopeptidase (oxytocinase)
    "CD157": 683,     # BST1/CD157 - involved in oxytocin secretion
})

# Example helper functions for working with the data structure
def get_all_genes():
    """Returns a set of all genes involved in neurotransmitter systems."""
    all_genes = set()
    
    def extract_genes(d):
        for k, v in d.items():
            if isinstance(v, list):
                all_genes.update(v)
            elif isinstance(v, dict):
                extract_genes(v)
    
    extract_genes(neurotransmitter_systems)
    return all_genes

def extract_genes(d):
    genes = set()
    for k, v in d.items():
        if isinstance(v, list):
            genes.update(v)
        elif isinstance(v, dict):
            extract_genes(v)
    return genes

def get_neurotransmitter_genes(neurotransmitter):
    """Returns a set of all genes associated with a specific neurotransmitter."""
    if neurotransmitter not in neurotransmitter_systems:
        return set()
    
    genes = set()
    def extract_genes(d):
        for k, v in d.items():
            if isinstance(v, list):
                genes.update(v)
            elif isinstance(v, dict):
                extract_genes(v)
    
    system = neurotransmitter_systems[neurotransmitter]
    extract_genes(system)
    return genes

def get_receptor_genes():
    """Returns a set of all receptor genes."""
    receptor_genes = set()
    
    for system in neurotransmitter_systems.values():
        if isinstance(system, dict) and "receptors" in system:
            def extract_receptor_genes(d):
                for k, v in d.items():
                    if isinstance(v, list):
                        receptor_genes.update(v)
                    elif isinstance(v, dict):
                        extract_receptor_genes(v)
            
            extract_receptor_genes(system["receptors"])
    
    return receptor_genes

def add_entrez_to_dict(d):
    new_dict = {}
    for k, v in d.items():
        if isinstance(v, list):
            new_dict[k] = [(gene, gene_to_entrez[gene]) for gene in v]
        elif isinstance(v, dict):
            new_dict[k] = add_entrez_to_dict(v)
        else:
            new_dict[k] = v
    return new_dict

# Create the new dictionary with Entrez IDs
neurotransmitter_systems_with_entrez = add_entrez_to_dict(neurotransmitter_systems)

# Helper functions for the new data structure
def get_all_genes_with_entrez():
    """Returns a dictionary of all genes with their Entrez IDs."""
    all_genes = {}
    
    def extract_genes(d):
        for k, v in d.items():
            if isinstance(v, list) and len(v) > 0 and isinstance(v[0], tuple):
                all_genes.update({gene: entrez for gene, entrez in v})
            elif isinstance(v, dict):
                extract_genes(v)
    
    extract_genes(neurotransmitter_systems_with_entrez)
    return all_genes

def get_neurotransmitter_genes_with_entrez(neurotransmitter):
    """Returns a dictionary of genes and their Entrez IDs for a specific neurotransmitter."""
    if neurotransmitter not in neurotransmitter_systems_with_entrez:
        return {}
    
    genes = {}
    def extract_genes(d):
        for k, v in d.items():
            if isinstance(v, list) and len(v) > 0 and isinstance(v[0], tuple):
                genes.update({gene: entrez for gene, entrez in v})
            elif isinstance(v, dict):
                extract_genes(v)
    
    system = neurotransmitter_systems_with_entrez[neurotransmitter]
    extract_genes(system)
    return genes

def get_neurotransmitter_gene_categories(neurotransmitter):
    """
    Returns two lists of genes for a given neurotransmitter:
    1. Genes involved in synthesis and transport
    2. Receptor genes
    
    Args:
        neurotransmitter (str): Name of the neurotransmitter system
        
    Returns:
        tuple: (synthesis_transport_genes, receptor_genes) where each is a list of gene symbols
    """
    if neurotransmitter not in neurotransmitter_systems:
        return [], []
        
    system = neurotransmitter_systems[neurotransmitter]
    
    synthesis_transport_genes = []
    receptor_genes = []
    
    # Extract synthesis and transport genes
    if 'synthesis_transport' in system:
        for category in system['synthesis_transport'].values():
            if isinstance(category, list):
                synthesis_transport_genes.extend(category)
            
    # Extract receptor genes
    if 'receptors' in system:
        for receptor_type in system['receptors'].values():
            if isinstance(receptor_type, list):
                receptor_genes.extend(receptor_type)
            elif isinstance(receptor_type, dict):
                for subtype in receptor_type.values():
                    receptor_genes.extend(subtype)
                    
    return synthesis_transport_genes, receptor_genes


def get_nt_set_bias(neurotransmitter_systems_dict, gene_to_entrez, HCT_Z2_MAT_HCT, Anno):
    genes = extract_genes(neurotransmitter_systems_dict)
    genes_entrez = [gene_to_entrez[gene] for gene in genes]
    genes_gw = dict(zip(genes_entrez, [1]*len(genes_entrez)))

    bias = AvgCTZ_Weighted(HCT_Z2_MAT_HCT, genes_gw, Method = 1)
    bias = AnnotateCTDat(bias, Anno)
    return bias

def get_nt_set_bias_v2(genes, gene_to_entrez, HCT_Z2_MAT_HCT, Anno):
    genes_entrez = [gene_to_entrez[gene] for gene in genes]
    genes_gw = dict(zip(genes_entrez, [1]*len(genes_entrez)))

    bias = AvgCTZ_Weighted(HCT_Z2_MAT_HCT, genes_gw, Method = 1)
    bias = AnnotateCTDat(bias, Anno)
    return bias

# Example 
serotonin_genes = get_neurotransmitter_genes_with_entrez("serotonin")
serotonin_genes