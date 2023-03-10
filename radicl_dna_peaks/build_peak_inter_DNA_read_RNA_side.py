#%%
import bioframe as bf
import numpy as np
import altair as alta
import subprocess
import pandas as pd
#%%
rna_radicl_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/raw/RNA/RADICL_iPSC_RNA.bed"
peak_read_tbl_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/peak_DNA_read_inter_tbl.tsv"

out_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/RADICL_peak_RNA.bed"

#%%
peak_read_inter_tbl = pd.read_csv(peak_read_tbl_file,delimiter="\t")

#%%
rna_df = pd.read_csv(rna_radicl_file,sep="\t",header=None)
#%%
peak_rna_inter_df=rna_df.loc[rna_df.iloc[:,3].isin(peak_read_inter_tbl.read_ID.to_list()),:]

#%%
peak_rna_inter_df.to_csv(out_file,sep="\t",header=False,index=False)
#%%
iter_csv = pd.read_csv(rna_radicl_file,sep="\t",header=None, iterator=True, chunksize=100_000)
df = pd.concat([chunk.loc[chunk.iloc[:,3].isin(peak_read_inter_tbl.read_ID.to_list()),:] for chunk in iter_csv])

# %%
