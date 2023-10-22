import yaml
import pandas as pd
from collections import OrderedDict
from pathlib import Path

run_config = '/home/zhangh5/work/Tapestri_project/cn-calling/train-normals/NB_train_normals.config.yaml'

with open(run_config, 'rb') as f:
    run_config = yaml.safe_load(f)
train_sample_rc_tsvs = OrderedDict(run_config['train_sample_rc_tsvs'])
output_dir = Path(run_config['output_dir'])
output_dir = output_dir / 'amp_ado'
if not output_dir.exists():
    output_dir.mkdir(parents=True, exist_ok=True)

# # per-sample
# for sample_i in train_sample_rc_tsvs.keys():
#     print(sample_i)
#     rc_df = pd.read_csv(train_sample_rc_tsvs[sample_i], sep='\t', index_col=0)
#     rc_df_binarized = rc_df.applymap(lambda x: 0 if x > 0 else 1)
#     amp_ado_stats = rc_df_binarized.sum(axis=0) / len(rc_df_binarized)

#     # amp_ado_stats.to_csv(output_dir / f'{sample_i}.amp_ado.csv', header=False)

# all 8-sample concatenated
rc_df = pd.concat(
    [pd.read_csv(train_sample_rc_tsvs[sample_i], sep='\t', index_col=0) for sample_i in train_sample_rc_tsvs.keys()], 
    axis=0
    )
rc_df_binarized = rc_df.applymap(lambda x: 0 if x > 0 else 1)
amp_ado_stats = rc_df_binarized.sum(axis=0) / rc_df_binarized.shape[0]

amp_ado_stats.to_csv(output_dir / f'8_normal_samples.amp_ado.csv', header=False)


