import os

configfile: "config/config.yaml"

# Pull variables from config
GENESETS = config["gene_sets"]
ANALYSIS = list(config["analysis_types"].keys())
MATCHING_CONFIGS = config["matching_configs"]
PROJDIR = config["ProjDIR"]
N_WEIGHTS = config["n_random_weights"]
N_PROCESSES = config.get("n_processes", 10)

# Function to auto-generate output directory based on matching parameters
def build_output_dir(config_name, sampling_mode, matched_params={}, set_level_params={}):
    """
    Auto-generate output directory name based on sampling mode and parameters.

    Examples:
        - random -> results/random
        - matched -> results/matched_CDS_WB_LOEUF
        - set_level_matched with propensity -> results/set_level_matched_CDS_WB_LOEUF_PropWeight_Tricubic
        - set_level_matched with best-of-N -> results/set_level_matched_CDS_WB_LOEUF_Best1000
    """
    base = "results"

    if sampling_mode == "random":
        return f"{base}/{config_name}"

    elif sampling_mode == "matched":
        # Gene-by-gene matching
        vars_str = "_".join(sorted(matched_params.get("matching_variables", ["CDS_length", "WB", "LOEUF"])))
        kernel = matched_params.get("kernel", "tricubic").capitalize()
        return f"{base}/{config_name}_{vars_str}_{kernel}"

    elif sampling_mode == "set_level_matched":
        # Set-level distribution matching
        vars_str = "_".join(sorted(set_level_params.get("matching_variables", ["CDS_length", "WB", "LOEUF"])))

        # Determine which algorithm is being used
        if set_level_params.get("use_propensity", False):
            kernel = set_level_params.get("propensity_kernel", "tricubic").capitalize()
            method = f"PropWeight_{kernel}"
        elif set_level_params.get("use_best_of_n", False):
            n_cand = set_level_params.get("n_candidates", 100)
            method = f"Best{n_cand}"
        elif set_level_params.get("use_sis", False):
            method = "SIS"
        else:
            method = "Rejection"

        return f"{base}/{config_name}_{vars_str}_{method}"

    else:
        # Fallback to config name
        return f"{base}/{config_name}"

# Build output directories for all matching configurations
OUTDIRS = {}
OUTDIR_TO_CONFIG = {}  # Reverse mapping: outdir -> config_name
for config_name, config_params in MATCHING_CONFIGS.items():
    sampling_mode = config_params.get("sampling_mode", "random")
    matched_params = config_params.get("matched_sampling", {})
    set_level_params = config_params.get("set_level_matched_sampling", {})

    outdir = build_output_dir(config_name, sampling_mode, matched_params, set_level_params)
    OUTDIRS[config_name] = outdir
    OUTDIR_TO_CONFIG[outdir] = config_name
    print(f"Config '{config_name}': {outdir}")

# Helper function to build full path to expression matrices
def get_expr_path(analysis):
    return os.path.join(PROJDIR, config["analysis_types"][analysis])

# Helper function to get config name from output path
def get_config_from_path(path):
    """Extract config name from output path by matching against OUTDIRS."""
    for outdir in OUTDIRS.values():
        if path.startswith(outdir):
            return OUTDIR_TO_CONFIG[outdir]
    raise ValueError(f"Could not determine config from path: {path}")

rule all:
    input:
        # Generate outputs for all combinations of config × analysis × geneset
        expand("{outdir}/{analysis}/{geneset}_bias_addP.csv",
               outdir=[OUTDIRS[c] for c in MATCHING_CONFIGS.keys()],
               analysis=ANALYSIS,
               geneset=GENESETS),
        expand("{outdir}/{analysis}/{geneset}_bias_addP_supercluster.csv",
               outdir=[OUTDIRS[c] for c in MATCHING_CONFIGS.keys()],
               analysis=ANALYSIS,
               geneset=GENESETS)



# Step 1: Generate random gene weights for each config × gene set × analysis type
rule generate_geneweights:
    input:
        expr=lambda wc: get_expr_path(wc.analysis),
        geneweights=lambda wc: config["gene_sets"][wc.geneset]["geneweights"]
    output:
        weights="{outdir}/{analysis}/null_weights/{geneset}_random_geneweights.csv"
    params:
        n=N_WEIGHTS,
        n_processes=N_PROCESSES,
        # Extract config-specific parameters dynamically
        sampling_mode=lambda wc: MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("sampling_mode", "random"),
        # Gene-by-gene matched sampling parameters
        kernel=lambda wc: MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("matched_sampling", {}).get("kernel", "tricubic"),
        bandwidth=lambda wc: MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("matched_sampling", {}).get("bandwidth", 10.0),
        matched_seed_arg=lambda wc: f"--seed {MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get('matched_sampling', {}).get('seed')}" if MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("matched_sampling", {}).get("seed") is not None else "",
        matched_vars=lambda wc: ",".join(MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("matched_sampling", {}).get("matching_variables", ["CDS_length", "WB", "LOEUF"])),
        # Set-level matched sampling parameters
        max_distance=lambda wc: MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("set_level_matched_sampling", {}).get("max_distance", 0.15),
        max_attempts=lambda wc: MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("set_level_matched_sampling", {}).get("max_attempts", 100),
        use_ks_test_arg=lambda wc: "--use_ks_test" if MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("set_level_matched_sampling", {}).get("use_ks_test", False) else "",
        ks_threshold=lambda wc: MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("set_level_matched_sampling", {}).get("ks_threshold", 0.3),
        set_level_seed_arg=lambda wc: f"--seed {MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get('set_level_matched_sampling', {}).get('seed')}" if MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("set_level_matched_sampling", {}).get("seed") is not None else "",
        set_level_vars=lambda wc: ",".join(MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("set_level_matched_sampling", {}).get("matching_variables", ["CDS_length", "WB", "LOEUF"])),
        # Propensity weighted parameters
        use_propensity_arg=lambda wc: "--use_propensity" if MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("set_level_matched_sampling", {}).get("use_propensity", False) else "",
        propensity_kernel=lambda wc: MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("set_level_matched_sampling", {}).get("propensity_kernel", "tricubic"),
        propensity_bandwidth=lambda wc: MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("set_level_matched_sampling", {}).get("propensity_bandwidth", 50.0),
        add_noise_arg=lambda wc: "--add_noise" if MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("set_level_matched_sampling", {}).get("add_noise", False) else "",
        noise_scale=lambda wc: MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("set_level_matched_sampling", {}).get("noise_scale", 5.0),
        relaxed_matching_arg=lambda wc: "--relaxed_matching" if MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("set_level_matched_sampling", {}).get("relaxed_matching", False) else "",
        loeuf_weight=lambda wc: MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("set_level_matched_sampling", {}).get("loeuf_weight", 0.5),
        # Best-of-N parameters
        use_best_of_n_arg=lambda wc: "--use_best_of_n" if MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("set_level_matched_sampling", {}).get("use_best_of_n", False) else "",
        n_candidates=lambda wc: MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("set_level_matched_sampling", {}).get("n_candidates", 100),
        # SIS (Sequential Importance Sampling) parameters
        use_sis_arg=lambda wc: "--use_sis" if MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("set_level_matched_sampling", {}).get("use_sis", False) else "",
        use_percentile_arg=lambda wc: "--use_percentile" if MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("set_level_matched_sampling", {}).get("use_percentile", False) else "",
        temperature=lambda wc: MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("set_level_matched_sampling", {}).get("temperature", 1.0),
        adaptive_temp_arg=lambda wc: "--adaptive_temp" if MATCHING_CONFIGS[get_config_from_path(wc.outdir)].get("set_level_matched_sampling", {}).get("adaptive_temp", False) else ""
    shell:
        """
        mkdir -p $(dirname {output.weights})
        if [ "{params.sampling_mode}" = "matched" ]; then
            python scripts/script_generate_geneweights.py \
                --WeightDF {input.geneweights} \
                --SpecMat {input.expr} \
                --n_sims {params.n} \
                --sampling_mode {params.sampling_mode} \
                --kernel {params.kernel} \
                --bandwidth {params.bandwidth} \
                --matching_variables {params.matched_vars} \
                --n_processes {params.n_processes} \
                {params.matched_seed_arg} \
                --outfile {output.weights}
        elif [ "{params.sampling_mode}" = "set_level_matched" ]; then
            python scripts/script_generate_geneweights.py \
                --WeightDF {input.geneweights} \
                --SpecMat {input.expr} \
                --n_sims {params.n} \
                --sampling_mode {params.sampling_mode} \
                --matching_variables {params.set_level_vars} \
                --max_distance {params.max_distance} \
                --max_attempts {params.max_attempts} \
                {params.use_ks_test_arg} \
                --ks_threshold {params.ks_threshold} \
                {params.use_propensity_arg} \
                --propensity_kernel {params.propensity_kernel} \
                --propensity_bandwidth {params.propensity_bandwidth} \
                {params.add_noise_arg} \
                --noise_scale {params.noise_scale} \
                {params.relaxed_matching_arg} \
                --loeuf_weight {params.loeuf_weight} \
                {params.use_best_of_n_arg} \
                --n_candidates {params.n_candidates} \
                {params.use_sis_arg} \
                {params.use_percentile_arg} \
                --temperature {params.temperature} \
                {params.adaptive_temp_arg} \
                --n_processes {params.n_processes} \
                {params.set_level_seed_arg} \
                --outfile {output.weights}
        else
            python scripts/script_generate_geneweights.py \
                --WeightDF {input.geneweights} \
                --SpecMat {input.expr} \
                --n_sims {params.n} \
                --sampling_mode {params.sampling_mode} \
                --n_processes {params.n_processes} \
                --outfile {output.weights}
        fi
        """

# Step 2: Run null bias calculation for each config × gene set × analysis type
rule compute_null_bias:
    input:
        expr=lambda wc: get_expr_path(wc.analysis),
        weights="{outdir}/{analysis}/null_weights/{geneset}_random_geneweights.csv"
    output:
        bias="{outdir}/{analysis}/null_bias/{geneset}_null_bias.csv"
    params:
        geneset="{geneset}"
    shell:
        """
        mkdir -p $(dirname {output.bias})
        python scripts/script_run_ctrl_sim.py \
            --SpecMat {input.expr} \
            --outfile {output.bias} \
            --mode human_ct_bias \
            --Ctrl_Genes_Fil {input.weights} \
        """

# Step 3: Run bias calculation for each config × gene set × analysis type
rule compute_bias_pvalue:
    input:
        expr=lambda wc: get_expr_path(wc.analysis),
        geneweights=lambda wc: config["gene_sets"][wc.geneset]["geneweights"],
        bias_null="{outdir}/{analysis}/null_bias/{geneset}_null_bias.csv"
    output:
        bias_final="{outdir}/{analysis}/{geneset}_bias_addP.csv",
        bias_supercluster="{outdir}/{analysis}/{geneset}_bias_addP_supercluster.csv"
    params:
        geneset="{geneset}"
    shell:
        """
        mkdir -p $(dirname {output.bias_final})
        python scripts/script_bias_cal.py \
            --SpecMat {input.expr} \
            --gw {input.geneweights} \
            --Bias_Null {input.bias_null} \
            --Bias_Out {output.bias_final} \
            --Bias_Out_Supercluster {output.bias_supercluster}
        """