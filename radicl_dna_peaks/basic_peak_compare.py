#%%
import pandas as pd
import numpy as np
import altair as alt
#%%
peak_file_a = "/home/vipink/Documents/FANTOM6/data/RADICL_data/IPSC/replicate1/results_beta/MACS/peaks/RADICL_DNA_tot_peaks.bed"
peak_file_b = "/home/vipink/Documents/FANTOM6/data/RADICL_data/Neuron/replicate1/results/MACS/peaks/RADICL_DNA_tot_peaks.bed"

#%%
peak_df_a = pd.read_csv(peak_file_a,header=None,delimiter="\t")
peak_df_a.columns = ['chrom','start','end','ID','score',"strand","enrichment","pvalue","qvalue","peak"]

peak_df_b = pd.read_csv(peak_file_b,header=None,delimiter="\t")
peak_df_b.columns = ['chrom','start','end','ID','score',"strand","enrichment","pvalue","qvalue","peak"]

#%%
peak_summary_df = (pd.concat([peak_df_a.assign(set='iPSC'),peak_df_b.assign(set='Neuron')])
                   .assign(peak_width=lambda df_:df_.end-df_.start))

# %%
alt.data_transformers.disable_max_rows()

dens_chart = (alt.Chart(peak_summary_df.assign(log_value=lambda df_:np.log10(df_.peak_width)))
.transform_density(
    'log_value',
    as_=['log_value', 'density'],
    groupby=['set']
).mark_line(size=1).encode(
    x=alt.X("log_value:Q",axis=alt.Axis(format='g', title='log10(peak width)')),
    y='density:Q',
    color='set:N'
))

# %%
dens_chart.save(f"/home/vipink/Documents/FANTOM6/litt/poster/img/peak_size_dens.html", scale_factor=1.0)

# %%
