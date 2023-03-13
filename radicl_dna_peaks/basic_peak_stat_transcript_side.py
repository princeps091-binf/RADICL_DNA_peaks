#%%
import pandas as pd
import bioframe as bf
import altair as alt
import numpy as np

#%%
rna_radicl_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/RADICL_peak_RNA.bed"
annotation_file = "/home/vipink/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed"

#%%
rna_df = pd.read_csv(rna_radicl_file,header=None,delimiter="\t")
rna_df.columns = ['chrom','start','end','read_ID','score','strand']

# %%
transcript_annotation_df = pd.read_csv(annotation_file,header=None,delimiter="\t")
#%%
transcript_annotation_df.columns = ['chrom','start','end',
                                    'ID','score','strand',
                                    'start_b','end_b','sym',
                                    'exon_count','exon_length','exon_start']
#%%
transcript_cluster_df = (bf.cluster(transcript_annotation_df,on=['strand'])
                         .loc[:,['chrom','cluster_start','cluster_end',
                                 'strand','cluster','ID']])
transcript_cluster_coord_df = (transcript_cluster_df
                               .groupby(['chrom','cluster_start','cluster_end','strand','cluster'])
                               .size()
                               .reset_index(name='transcript_count')
                               .rename(columns ={'cluster_start': 'start',
                                                 'cluster_end':'end'}))
# %%
(bf.count_overlaps(rna_df,transcript_cluster_coord_df)
 .query('count<1'))
# %%
(bf.count_overlaps(transcript_cluster_coord_df,rna_df)
 .query('count>0')
 .sort_values('count'))

# %%
