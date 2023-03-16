#%%
import pandas as pd
import bioframe as bf
import numpy as np
from joblib import Parallel, delayed
import multiprocessing
#%%
intron_out_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/intron_coord.bed"
exon_out_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/exon_coord.bed"

#%%
transcript_annotation_file="/home/vipink/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed"

transcript_annotation_df = pd.read_csv(transcript_annotation_file,
                                       header=None,delimiter="\t")

transcript_annotation_df.columns = ['chrom','start','end',
                                    'ID','score','strand',
                                    'start_b','end_b','sym',
                                    'exon_count','exon_length','exon_start']

#%%
transcript_annotation_df = (transcript_annotation_df
 .assign(length_array=lambda df_ : df_.exon_length
         .str
         .split(',')
         .apply(lambda df_ : np.array(df_,dtype=np.int64)))
 .assign(start_pos_array=lambda df_ : df_.exon_start
         .str
         .split(',')
         .apply(lambda x : np.array(x,dtype=np.int64))))
#%%
exon_df = (transcript_annotation_df
 .assign(exon_start_array= lambda df_ : df_.start + df_.start_pos_array)
 .assign(exon_end_array= lambda df_ : df_.length_array + df_.exon_start_array)
 .loc[:,('chrom','start','end','strand','ID','exon_start_array','exon_end_array')]
 .explode(['exon_start_array', 'exon_end_array']))
#%%
tmp_set = ["CATG00000114975.1|MICT00000384061.1",'ENSG00000225880.4|FTMT20100027365.1']
transcript_test = (exon_df
                   .query(" @tmp_set == ID"))
#%%
def produce_intron(df_,name):
    tmp_transcript_df = df_.loc[:,('chrom','start','end','strand')].drop_duplicates()
    tmp_exon_df = df_.loc[:,('chrom','strand','exon_start_array','exon_end_array')]
    tmp_exon_df = (tmp_exon_df
                   .assign(start= lambda df_: df_.loc[:,'exon_start_array'].astype(int),
                           end=lambda df_: df_.loc[:,'exon_end_array'].astype(int))
                   .loc[:,('chrom','start','end','strand')])
    return bf.subtract(tmp_transcript_df,tmp_exon_df).assign(ID=name)

def applyParallel(dfGrouped, func):
    retLst = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(func)(group,name) for name, group in dfGrouped)
    return pd.concat(retLst)

#%%
intron_df = (exon_df
 .groupby("ID")
 .apply(lambda df_:produce_intron(df_))
 .reset_index())

# %%
intron_df = applyParallel(exon_df.groupby("ID"), produce_intron)

intron_df.to_csv(intron_out_file,
                 sep="\t",
                 header=True,
                 index=False)

#%%
(exon_df
 .loc[:,['chrom','exon_start_array','exon_end_array','strand','ID']]
 .rename(columns={'exon_start_array':'start',
                  'exon_end_array':'end'})
 .to_csv(exon_out_file,
         sep="\t",
         header=True,
         index=False))
# %%
