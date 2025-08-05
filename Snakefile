import os

configfile: "config/config.yaml"

# Pull variables from config
GENESETS = config["gene_sets"]
ANALYSIS = list(config["analysis_types"].keys())
PROJDIR = config["ProjDIR"]
OUTDIR = config["output_dir"]
N_WEIGHTS = config["n_random_weights"]

# Helper function to build full path to expression matrices
def get_expr_path(analysis):
    return os.path.join(PROJDIR, config["analysis_types"][analysis])

rule all:
    input:
        expand(f"{OUTDIR}/{{analysis}}/{{geneset}}_bias_addP.csv",
               analysis=ANALYSIS,
               geneset=GENESETS)



# Step 1: Generate random gene weights for each gene set × analysis type
rule generate_geneweights:
    input:
        expr=lambda wc: get_expr_path(wc.analysis),
        geneweights=lambda wc: config["gene_sets"][wc.geneset]["geneweights"]
    output:
        weights=f"{OUTDIR}/{{analysis}}/null_weights/{{geneset}}_random_geneweights.csv"
    params:
        n=N_WEIGHTS
    shell:
        """
        mkdir -p $(dirname {output.weights})
        python scripts/script_generate_geneweights.py \
            --WeightDF {input.geneweights} \
            --SpecMat {input.expr} \
            --n_sims {params.n} \
            --outfile {output.weights}
        """

# Step 2: Run null bias calculation for each gene set × analysis type
rule compute_null_bias:
    input:
        expr=lambda wc: get_expr_path(wc.analysis),
        weights=f"{OUTDIR}/{{analysis}}/null_weights/{{geneset}}_random_geneweights.csv"
    output:
        bias=f"{OUTDIR}/{{analysis}}/null_bias/{{geneset}}_null_bias.csv"
    params:
        geneset="{geneset}"
    shell:
        """
        mkdir -p $(dirname {output.bias})
        python scripts/script_run_ctrl_sim.py \
            --SpecMat {input.expr} \
            --WeightDF {input.weights} \
            --outfile {output.bias} \
            --mode human_ct_bias \
            --Ctrl_Genes_Fil {input.weights} \
        """

# Step 3: Run bias calculation for each gene set × analysis type
rule compute_bias_pvalue:
    input:
        expr=lambda wc: get_expr_path(wc.analysis),
        geneweights=lambda wc: config["gene_sets"][wc.geneset]["geneweights"],
        bias_null=f"{OUTDIR}/{{analysis}}/null_bias/{{geneset}}_null_bias.csv"
    output:
        bias_final=f"{OUTDIR}/{{analysis}}/{{geneset}}_bias_addP.csv"
    params:
        geneset="{geneset}"
    shell:
        """
        mkdir -p $(dirname {output.bias_final})
        python scripts/script_bias_cal.py \
            --SpecMat {input.expr} \
            --gw {input.geneweights} \
            --Bias_Null {input.bias_null} \
            --Bias_Out {output.bias_final}
        """