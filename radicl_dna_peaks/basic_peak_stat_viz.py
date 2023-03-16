#%%
import pandas as pd
import altair as alt
import numpy as np

#%%
peak_inter_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/peak_DNA_read_inter_tbl.tsv"
peak_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA/MACS/peaks/RADICL_DNA_tot_peaks.bed"

#%%
# Get peak to reads
peak_read_inter_tbl = pd.read_csv(peak_inter_file,delimiter="\t")
#%%
peak_df = pd.read_csv(peak_file,header=None,delimiter="\t")
peak_df.columns = ['chrom','start','end','ID','score',"strand","enrichment","pvalue","qvalue","peak"]

# %%
peak_summary_df = (peak_read_inter_tbl
 .groupby(['chrom','start','end'])
 .size()
 .reset_index(name='read_count')
 .assign(read_density=lambda df_: df_.read_count/(df_.end-df_.start),
         peak_width=lambda df_:df_.end-df_.start))
# %%
alt.data_transformers.disable_max_rows()

peak_summary_df = (peak_summary_df
                   .assign(log_read=lambda df_:np.log10(df_.read_density),
                           log_width=lambda df_:np.log10(df_.peak_width)))
(alt.Chart(peak_summary_df)
.transform_density(
    'log_read',
    as_=['log_read', 'density'])
.mark_line()
.encode(
    alt.X("log_read:Q"),
    y='density:Q'))
# %%
(alt.Chart(peak_summary_df)
.mark_point(
    size=0.1,
    filled=True,
    opacity=0.5
)
.encode(
    alt.X("log_width:Q"),
    alt.Y('log_read:Q')))

# %%
