#%%
import pandas as pd
import bioframe as bf
from joblib import Parallel, delayed
import multiprocessing

#%%
exon_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/exon_coord.bed"
intron_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/intron_coord.bed"
RNA_read_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/raw/RNA/RADICL_iPSC_RNA.bed"
#%%
exon_df = pd.read_csv(exon_file,
                      delimiter="\t")
rna_df = pd.read_csv(RNA_read_file,header=None,delimiter="\t")

# %%
rna_df.columns = ['chrom','start','end','ID','score','strand']
# %%
exon_df = bf.count_overlaps(exon_df,rna_df,on=['strand'])
# %%
transcript_df = (exon_df
 .groupby("ID")
 .agg(chrom=('chrom','first'),
      start =('start', min ), 
      end =('end', max ),
      strand=('strand','first'))
 .reset_index())

# %%
def applyParallel(dfGrouped, func):
    retLst = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(func)(group,name) for name, group in dfGrouped)
    return pd.concat(retLst)
def read_count_fn(df_):
    bf.count_overlaps(df_,rna_df.query(f"chrom == '{df_.name}'"),on=['strand'])

# %%
(transcript_df
 .groupby('chrom')
 .apply(lambda df_:read_count_fn(df_)))
# %%
