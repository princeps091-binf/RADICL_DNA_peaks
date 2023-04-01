#%%
import pandas as pd
import altair as alt
import numpy as np
import bioframe as bf
#%%
rep1_peak_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/IPSC/replicate1/results_beta/MACS/peaks/RADICL_DNA_tot_peaks.bed"
rep2_peak_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/IPSC/replicate2/results/MACS/peaks/RADICL_DNA_tot_peaks.bed"

# %%
rep1_peak_df = pd.read_csv(rep1_peak_file,header=None,delimiter="\t")
rep1_peak_df.columns = ['chrom','start','end','ID','score',"strand","enrichment","pvalue","qvalue","peak"]

rep2_peak_df = pd.read_csv(rep2_peak_file,header=None,delimiter="\t")
rep2_peak_df.columns = ['chrom','start','end','ID','score',"strand","enrichment","pvalue","qvalue","peak"]

# %%
(bf.coverage(rep1_peak_df,rep2_peak_df)
.agg({'coverage':['sum']}))
# %%
(bf.coverage(rep2_peak_df,rep1_peak_df)
 .assign(w=lambda df_: df_.end - df_.start)
 .agg({'coverage':['sum'],'w':['sum']})
 .assign(r=lambda df_:df_.coverage/df_.w))
# %%
(bf.merge(pd.concat([rep2_peak_df,rep1_peak_df]))
 .assign(w=lambda df_: df_.end - df_.start)
 .agg({'w':['sum']}))
# %%
