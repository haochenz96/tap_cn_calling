# Copy number calling for Tapestri data

@Palash Sashittal, @Haochen Zhang

Use negative-binomial (NB) clustering to assign integer copy number (CN) state for each gene, for each single cell in Tapestri data. Further, amplicon-level homozygous-deletion (ploidy=0) can be detected.

## Install
It is recommended to use Conda/Mamba to install Snakemake and other dependencies. A script is provided in `setup_smk_env.sh`

## Run the program
The program runs in two phases:

### Phase 1: gene-level CN calling without homdel
    

Please see `examples/M04-sc_cn_calling_config.yaml` for details on configuring inputs/parameters.

- To run:
    ```
    snakemake \
        -s cn_call_no_homdel.smk \
        --configfile examples/M04-sc_cn_calling_config.yaml \
        --profile lsf \
        --conda-prefix /home/zhangh5/work/Tapestri_batch2/analysis/conda \
        --quiet rules
    ```
    
- *Output:*

    Three folders will be generated in the working directory:
    ```
    cn_call_no_homdel/intermediate_results/
    cn_call_no_homdel/solutions/
    cn_call_no_homdel/outputs/
    ```

    `intemediate_results/` contains EM run output for each `nclones` and each random seed.

    `solutions/` contains the best (lowest Bayesian Information Criterion) solution for each `nclones`

    `outputs/` contains a lineplot of BIC across different `nclones` and symlink to the solution with the best `nclones`.
    
### Phase 2: gene-level CN calling with amplicon-level homdel detection
    
Importantly, the `best_nclones` parameter needs to be determined from phase1's BIC plot and set in the configuration file.

- To run:
    ```
    snakemake \
        -s cn_call_with_homdel.smk \
        --configfile examples/M04-sc_cn_calling_config.yaml \
        --profile lsf \
        --conda-prefix /home/zhangh5/work/Tapestri_batch2/analysis/conda \
        --quiet rules
    ```

- *Output:*

    Three folders will be generated in the working directory:
    ```
    cn_call_with_homdel/intermediate_results/
    cn_call_with_homdel/solutions/
    cn_call_with_homdel/outputs/
    ```

    `intemediate_results/` contains EM run output for each `best_nclones`, if multiple were set.

    `solutions/` contains the compiled solution files for each `best_nclones`. These include:
    ```
    {cohort_name}-homdel-nclones={best_nclones}_start_from_best_sol_solution.amp_clone_profiles.csv'
    
    {cohort_name}-homdel-nclones={best_nclones}_start_from_best_sol_solution.cell_assignments.csv'
    ```
    `outputs/` contains:
    - a plot of each CN clone's per-amplicon ploidy state.
    - a plot of each sample's CN clone composition.

## Run Parameters

- panel_nseeds: number of different random seeds to try for EM with a given `nclones`
- panel_maxcn: maximum copy number state allowed 
- panel_amplicon_parameters: use `train-normals/train-combined_8_normals/NB_train-combined_8_normals-results.gene_added.csv`
- nclones (to try): a list of positive integers to use as the fixed number of clusters for each EM run.