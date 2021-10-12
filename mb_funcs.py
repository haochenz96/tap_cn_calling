import pandas as pd
import numpy as np
import loompy
import matplotlib.pyplot as plt

class analyzer():
    
    def __init__(self, ds):
        self.ds = ds
        self.generate_position_dataframe()
    
    def generate_position_dataframe(self):
        df_pos = pd.DataFrame({'pos': list(self.ds.ra['POS']),
                               'chrom': list(self.ds.ra['CHROM']),
                               'amplicon': list(self.ds.ra['amplicon']),
                               'ref': list(self.ds.ra['REF']),
                               'alt': list(self.ds.ra['ALT'])})

        df_pos['ref_len'] = df_pos['ref'].apply(len)
        df_pos['alt_len'] = df_pos['alt'].apply(len)
        df_pos['normal'] = (df_pos['ref_len'] == 1) & (df_pos['alt_len'] == 1)
        df_pos['index'] = df_pos.index

        self.df_unique_pos = df_pos.groupby(['pos', 'amplicon']).first().reset_index().sort_values('index')
        self.amplicon_list = list(self.df_unique_pos['amplicon'].unique())

    def plot_amplicon_coverage(self, curr_amplicon, nprobe = 5, trim_perc = 0, read_threshold = 0):
        amplicon_idx = self.amplicon_list.index(curr_amplicon)
        selected_indices = self.df_unique_pos[self.df_unique_pos['amplicon'] == curr_amplicon]['index'].values

        npos = len(selected_indices)
        if trim_perc > 0:
            cutoff = int(npos // (100 / trim_perc))
        else:
            cutoff = 0

        a = self.ds.layers['DP'][selected_indices[cutoff:npos-cutoff], :]    

        median_depths = np.median(a, axis = 0)
        a_filtered = a[:, median_depths > read_threshold]
        a_filtered_normalized = a_filtered / np.median(a_filtered, axis = 0)
        print(f"rejected {a.shape[1] - a_filtered_normalized.shape[1]} cells for amplicon {curr_amplicon}")

        np.random.seed(int(curr_amplicon.lstrip('AMPL')))

        fig, ax = plt.subplots(1,1,figsize=(10,4))

        if nprobe > 0:
            cell_idx = np.random.randint(a_filtered_normalized.shape[1], size = nprobe)    
            plt.plot(a_filtered_normalized[:,cell_idx], linewidth=3)
        else:
            cell_idx = np.random.randint(a_filtered_normalized.shape[1], size = nprobe)    
            plt.plot(a_filtered_normalized[:, :], 'k', linewidth=3, alpha = 0.05)        
        plt.title(f'{curr_amplicon} coverage for {nprobe} random cells from organoid', fontsize=20)
        plt.gca().set_ylabel('normalized depth')
        plt.gca().set_xlabel('index of position in the amplicon')
        
def get_position_dataframe_from_loom(ds):
    df_pos = pd.DataFrame({'pos': list(ds.ra['POS']),
                           'chrom': list(ds.ra['CHROM']),
                           'amplicon': list(ds.ra['amplicon']),
                           'ref': list(ds.ra['REF']),
                           'alt': list(ds.ra['ALT'])})

    df_pos['ref_len'] = df_pos['ref'].apply(len)
    df_pos['alt_len'] = df_pos['alt'].apply(len)
    df_pos['normal'] = (df_pos['ref_len'] == 1) & (df_pos['alt_len'] == 1)
    df_pos['index'] = df_pos.index

    return df_pos.groupby(['pos', 'amplicon']).first().reset_index().sort_values('index')