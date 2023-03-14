# set up Snakemake environment

# # ----- install mambaforge if needed -----
# wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
# bash Mambaforge-$(uname)-$(uname -m).sh

# create Snakemake environment
mamba create -c conda-forge -c bioconda -n snakemake snakemake

# install plotly and kaleido for plotting
mamba install -n snakemake -c plotly plotly python-kaleido