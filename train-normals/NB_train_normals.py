# %%
import sys, argparse
from pathlib import Path
import yaml
from collections import OrderedDict
import statsmodels.api as sm
import pandas as pd
import numpy as np
from IPython import embed
import pickle

# %% Define functions
def fit_NB_model(train_array, X, exposure, maxiter = 100):
    '''
    args:
    train_array: type = np.ndarray
        needs to be an 1-D array
    
    X: type = pd.DataFrame

    exposure: type = pd.Series

    '''
    # y = df_tsv[amplicon]
    nb_model = sm.NegativeBinomial(
        train_array, 
        X, 
        exposure = exposure
        )
    # nb_model = sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=1.0/fixed_dispersion))
    nb_res = nb_model.fit(disp=False, maxiter = maxiter)
    return nb_res
# fit_NB_model = np.vectorize(fit_NB_model, excluded = ['X', 'exposure', 'maxiter'])

def extract_NB_params(nb_res, amplicon_name, amplicon_rc, ):
    '''
    Args:
        nb_res_array: scifit results object

        amplicon_name: str
            the amplicon ID (e.g. 'AMPL41099')

        amplicon_rc: type = np.array
            needs to be an 1-D array holding the read counts for this particular amplicon

    Returns: type = list
        [amplicon_name, nb_res.mle_retvals['converged'], amplicon_factor, alpha, beta_zero, beta_one, 'pearson', residual_mean, residual_var]

    '''
    
    # amplicon_rc = df_tsv[amplicon]
    
    beta_zero = nb_res.params['constant']
    # beta_one = nb_res.params['total']
    beta_one = 1.0
    alpha = nb_res.params['alpha']
    # alpha = 1 / fixed_dispersion
    nb_mu = np.exp(beta_zero) * (amplicon_rc**beta_one)
    nb_var = nb_mu + nb_mu**2 * alpha
    amplicon_factor = np.exp(beta_zero) / 2
    pearson_residual = (amplicon_rc - nb_mu) / np.sqrt(nb_var)
    residual_mean = np.mean(pearson_residual)
    residual_var = np.var(pearson_residual)
    return np.array(
        [amplicon_name, nb_res.mle_retvals['converged'], amplicon_factor, alpha, beta_zero, beta_one, 'pearson', residual_mean, residual_var]
    )

# df_nb_residual = pd.DataFrame(
#     nb_residual_data, 
#     columns = ['amplicon', 'amplicon_factor', 'alpha', 'beta_zero', 'beta_one', 'method', 'mean', 'variance']
#     )

# %% Main
# def main(args):
# ----- read in inputs -----
# with open(args.run_config, 'rb') as f:
#     run_config = yaml.safe_load(f)
run_config_f = "/juno/work/iacobuzc/haochen/Tapestri_project/cn-calling/train-normals/NB_train_normals.config.yaml"
with open(run_config_f, 'rb') as f:
    run_config = yaml.safe_load(f)
train_sample_rc_tsvs = OrderedDict(run_config['train_sample_rc_tsvs'])
output_dir = Path(run_config['output_dir'])
if not output_dir.exists():
    output_dir.mkdir(parents=True, exist_ok=True)
run_topic = run_config['run_topic']
amp_gene_map_f = Path(run_config['amp_gene_map_f'])

# if multiple samples are provided, combine them into one dataframe
if len(train_sample_rc_tsvs) > 1:
    combined_df = pd.concat(
        [pd.read_csv(rc_i, sep='\t', index_col=0) for rc_i in train_sample_rc_tsvs.values()], 
        keys=[i for i in train_sample_rc_tsvs.keys()],
        names = ['sample_name', 'cell_name']
        )
    merged_index = [ f'{i}-{j}' for i, j in combined_df.index ]
    combined_df.index = merged_index
else:
    combined_df = pd.read_csv(list(train_sample_rc_tsvs.values())[0], sep='\t', index_col=0)
    # combined_df.index = [ f'{i}-{j}' for i, j in combined_df.index ]

sc_total_rc = combined_df.sum(axis = 1)
# X = np.log(sc_total_rc).to_frame(name='total')
X = pd.DataFrame(
    [1.0] * len(sc_total_rc), 
    index = sc_total_rc.index, 
    columns=['constant']
    )

# %% Fit model
# ----- fit NB model -----
nb_res_results = pd.Series(
    [fit_NB_model(combined_df[amplicon_i], X, sc_total_rc) for amplicon_i in combined_df.columns],
    index = combined_df.columns
    ) # N x 1 array storing the training results of each amplicon
print('[SUCCESS] NB training finished successfully')

# %%
# ----- pickle the NB results -----
output_f = Path(output_dir) / f'{run_topic}-results.pkl'
with open(output_f, 'wb') as output:
    pickle.dump(nb_res_results, output)
    print(f'[INFO] NB training results pickled to {output_f}')

nb_res_results_df = pd.DataFrame(
    index = nb_res_results.index,
    columns = ['amplicon', 'converged', 'amplicon_factor', 'alpha', 'beta_zero', 'beta_one', 'method', 'mean', 'variance']
    )
for amplicon_i in nb_res_results.index:
    nb_res_results_df.loc[amplicon_i, :] = extract_NB_params(nb_res_results[amplicon_i], amplicon_i, combined_df[amplicon_i])
    if nb_res_results_df.loc[amplicon_i, 'converged'] == False:
        print('=' * 8)
        print(f'[WARNING] {amplicon_i} did not converge')

ref_df = pd.read_csv(amp_gene_map_f, sep='\t')
ref_df.index = ref_df['amplicon_number']
nb_res_results_df['gene'] = ref_df.loc[
    nb_res_results_df['amplicon'],
    'gene_name'
] # get mapped gene names
nb_res_results_df.to_csv(Path(output_dir) / f'{run_topic}-results.csv', index=False, header=True)

print('[INFO] NB training results saved to:', Path(output_dir) / f'{run_topic}-results.tsv')
# embed()

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--run_config', type=str, help='input config yaml file')
#     # parser.add_argument('--output_dir', type=str, help='output directory')

#     args = parser.parse_args(None if sys.argv[1:] else ['-h'])

#     main(args)
# %%
