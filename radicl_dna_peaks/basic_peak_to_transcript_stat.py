#%%
import pandas as pd
import bioframe as bf
import altair as alt
import numpy as np
#%%
peak_to_transcript_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/python_project_result/DNA_peak_analysis/peak_to_transcript_tbl.tsv"
black_list_file = "/home/vipink/Documents/FANTOM6/data/annotation/hg38-blacklist.v2.bed"

#%%
peak_to_transcript_df=pd.read_csv(peak_to_transcript_file,sep='\t')
black_list_df=pd.read_csv(black_list_file,sep="\t",header=None)
black_list_df.columns = ['chrom','start','end','label']
#%%
#filter out black-list peaks
clean_peak_to_transcript_df = (bf.count_overlaps(peak_to_transcript_df
 .rename(columns={'chrom_peak':'chrom','start_peak':'start','end_peak':'end'}),black_list_df)
 .query('count < 1'))
# %%
peak_summary_df = (clean_peak_to_transcript_df
 .assign(intra=lambda df_:df_.chrom == df_.chrom_transcript)
 .groupby(['chrom','start','end'])
 .agg(intra_prop=('intra','mean'),
      read_count=('read_ID','count'))
 .reset_index()
 .assign(width=lambda df_:df_.end-df_.start)
 .sort_values('intra_prop'))

#%%
alt.data_transformers.disable_max_rows()

(alt.Chart(peak_summary_df)
.transform_density(
    'intra_prop',
    as_=['intra_prop', 'density'])
.mark_line()
.encode(
    alt.X("intra_prop:Q"),
    y='density:Q'))

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
peak_transcript_summary_df = (clean_peak_to_transcript_df
 .assign(intra=lambda df_:df_.chrom == df_.chrom_transcript)
 .groupby(['chrom','start','end'])
 .agg(intra_prop=('intra','mean'),
      gene_io=('genic',lambda df_: df_.isin(['genic']).mean()),
      read_count=('read_ID','count'))
 .reset_index()
 .assign(width=lambda df_:df_.end-df_.start)
 .assign(read_dens=lambda df_: df_.read_count/df_.width)
 .sort_values('read_count'))

#%%
(alt.Chart(peak_transcript_summary_df.assign(log_size=lambda df_:np.log10(df_.read_dens)))
.transform_density(
    'intra_prop',
    as_=['intra_prop', 'density'])
.mark_line()
.encode(
    alt.X("intra_prop:Q"),
    y='density:Q'))


#%%
(alt.Chart(peak_transcript_summary_df.assign(log_size=lambda df_:np.log10(df_.width)))
.mark_point(size=0.1,opacity=0.1)
.encode(
    alt.X("gene_io:Q"),
    alt.Y('log_size:Q')))

# %%
peak_transcript_summary_df = (clean_peak_to_transcript_df
                              .assign(target_ID=lambda df_:df_.ID
                                      .mask(df_.ID.isnull(),'intergenic'))
                              .groupby(['chrom','start','end','target_ID'])
                              .agg(read_count=('read_ID','count'))
                              .reset_index()
                              .sort_values('read_count')
                              .assign(width=lambda df_:df_.end-df_.start))

#%%
peak_summary = (clean_peak_to_transcript_df
 .loc[:,['chrom','start','end','read_ID','genic']]
 .drop_duplicates()
 .groupby(['chrom','start','end'])
 .agg(gene_io=('genic',lambda df_: df_.isin(['genic']).mean()),
      peak_nread=('read_ID','size'))
 .reset_index())

#%%
peak_transcript_summary_df = (peak_transcript_summary_df
 .merge(peak_summary)
 .assign(read_prop=lambda df_:df_.read_count/df_.peak_nread))

#%%
max_transcript_cov_df = (peak_transcript_summary_df.loc[(peak_transcript_summary_df
 .query("target_ID != 'intergenic'")
 .groupby(['chrom','start','end'])
 .read_prop
 .idxmax())]
 .sort_values('read_count'))

#%%
alt.data_transformers.disable_max_rows()

(alt.Chart(max_transcript_cov_df.assign(log_value=lambda df_:np.log10(df_.read_count)))
.transform_density(
    'log_value',
    as_=['log_value', 'density'],
).mark_line(size=1).encode(
    x=alt.X("log_value:Q",axis=alt.Axis(format='g', title='log10(peak width)')),
    y='density:Q',
))

#%%
(peak_transcript_summary_df
 .merge(peak_summary)
 .assign(read_prop=lambda df_:df_.read_count/df_.peak_nread)
 .groupby('target_ID')
 .agg(min=('read_count',min),
      max=('read_count',max),
      mean=('read_count','median'),
      minp=('read_prop',min),
      maxp=('read_prop',max),
      meanp=('read_prop','mean'),
      npeak=('chrom','count'))
  .query("npeak > 1")
  .sort_values('meanp'))
# %%
(peak_transcript_summary_df
 .groupby('target_ID')
 .agg(min=('read_count',min),
      max=('read_count',max),
      mean=('read_count','median'),
      npeak=('chrom_peak','count'))
 .query("npeak < 1_000")
 .sort_values('npeak'))
# %%
