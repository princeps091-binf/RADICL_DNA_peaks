#%%
import bioframe as bf
import numpy as np
import altair as alta
import subprocess
import pandas as pd

#%%
peak_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA/MACS/peaks/RADICL_DNA_tot_peaks.bed"
DNA_read_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/raw/DNA/RADICL_iPSC_DNA.bed"
peak_read_inter_tbl_out_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/peak_DNA_read_inter_tbl.tsv"
peak_df = pd.read_csv(peak_file,header=None,delimiter="\t")
DNA_reads_tbl = pd.read_csv(DNA_read_file,header=None,delimiter="\t")

DNA_reads_tbl.columns = ['chrom','start','end','ID','score']
peak_df.columns = ['chrom','start','end','ID','score',"char","score_a","score_b","score_c","score_d"]


# %%
peak_read_inter_idx = bf.overlap(peak_df, DNA_reads_tbl, how='inner', suffixes=('_1','_2'),return_index=True,return_input=False)

# %%
 
peak_read_inter_tbl = (peak_df.loc[peak_read_inter_idx.index_1,["chrom","start","end"]].reset_index(drop=True)
 .assign(read_ID=DNA_reads_tbl.loc[peak_read_inter_idx.index_2,["ID"]].reset_index(drop=True)))

# %%
peak_read_inter_tbl.to_csv(peak_read_inter_tbl_out_file,
                           sep="\t",
                           header=True,
                           index=False)