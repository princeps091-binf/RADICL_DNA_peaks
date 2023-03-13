#%%
import pandas as pd
import bioframe as bf
import numpy as np

#%%
annotation_file = "/home/vipink/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed"
out_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/transcript_cluster.bed"
#%%
transcript_annotation_df = pd.read_csv(annotation_file,header=None,delimiter="\t")

transcript_annotation_df.columns = ['chrom','start','end',
                                    'ID','score','strand',
                                    'start_b','end_b','sym',
                                    'exon_count','exon_length','exon_start']


# %%
(bf.cluster(transcript_annotation_df,on=['strand'])
 .loc[:,['chrom','cluster_start','cluster_end','strand','cluster','ID']]
 .groupby('cluster')
 .size()
 .reset_index(name='counts')
 .sort_values('counts')
 .assign(log_value=lambda df_: np.log10(df_.loc[:,['counts']]))
 .loc[:,'log_value']
 .plot.density ()
)
# %%
transcript_cluster_df = (bf.cluster(transcript_annotation_df,on=['strand'])
                         .loc[:,['chrom','cluster_start','cluster_end',
                                 'strand','cluster','ID']])
transcript_cluster_coord_df = (transcript_cluster_df
                               .rename(columns={'cluster_start':'start',
                                        'cluster_end':'end'})
                               .groupby(['chrom','start','end','strand','cluster'])
                               .size()
                               .reset_index(name="transcript_count"))
# %%
transcript_cluster_coord_df.to_csv(out_file,
                           sep="\t",
                           header=True,
                           index=False)

# %%
