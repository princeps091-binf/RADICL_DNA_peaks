#%%
import pandas as pd
import bioframe as bf
import altair as alt
import numpy as np
#%%
peak_to_transcript_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/peak_to_transcript_tbl.tsv"

#%%
peak_to_transcript_df=pd.read_csv(peak_to_transcript_file,sep='\t')
# %%
peak_summary_df = (peak_to_transcript_df
 .assign(intra=lambda df_:df_.chrom_peak == df_.chrom_transcript)
 .groupby(['chrom_peak','start_peak','end_peak'])
 .agg(intra_prop=('intra','mean'),
      read_count=('read_ID','count'))
 .reset_index()
 .assign(width=lambda df_:df_.end_peak-df_.start_peak)
 .sort_values('intra_prop'))
# %%
alt.data_transformers.disable_max_rows()

(alt.Chart(peak_summary_df.assign(log_size=lambda df_:np.log10(df_.width)))
.transform_density(
    'log_size',
    as_=['log_size', 'density'])
.mark_line()
.encode(
    alt.X("log_size:Q"),
    y='density:Q'))

# %%
(alt.Chart(peak_summary_df.assign(log_size=lambda df_:np.log10(df_.read_count)))
.mark_point(size=0.1,opacity=0.1)
.encode(
    alt.X("intra_prop:Q"),
    alt.Y('log_size:Q')))

# %%
peak_transcript_summary_df = (peak_to_transcript_df
 .assign(intra=lambda df_:df_.chrom_peak == df_.chrom_transcript)
 .groupby(['chrom_peak','start_peak','end_peak','ID'])
 .agg(read_count=('read_ID','count'))
 .reset_index()
 .assign(width=lambda df_:df_.end_peak-df_.start_peak)
 .sort_values('read_count'))

# %%
