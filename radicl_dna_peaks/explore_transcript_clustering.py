#%%
import pandas as pd
import bioframe as bf
import numpy as np

#%%
annotation_file = "/home/vipink/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed"

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
