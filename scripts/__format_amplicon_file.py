# %%
import pandas as pd

panel_f = "/data/iacobuzc/haochen/Tapestri_project/panels/Tapestri-Designer-results-Myeloid/Myeloid.exon_coverage.bed"
panel = pd.read_csv(panel_f, sep="\t", header=None)
coi = [0,1,2,3,7]
panel = panel[coi]
col_names = ["chrom","min_pos","max_pos","amplicon","gene"]
panel.columns = col_names
panel.to_csv("/data/iacobuzc/haochen/Tapestri_Rocio/Myeloid_panel.amplicon_coords.csv", header=True, index=False)



# %%
