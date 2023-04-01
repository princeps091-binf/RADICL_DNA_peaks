#%%
import pandas as pd
import os
import bioframe as bf
import re
#iPSC
#47843486 total reads in RADICL
#33373460 reads intersecting transcript
#5636408 reads intersecting enhancers
#autosome
#46041920  total reads in RADICL
#32517574 reads intersecting transcript
#5556254 reads intersecting enhancers

#Neuron
#74004761 total reads
#53798140 reads intersecting transcript
#9601693 reads intersecting enhancers
trx_read_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/DNA_read_transcript_intersect.bed"
enh_read_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/DNA_read_enhancer_intersect.bed"
annotation_file = "/home/vipink/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed"
enh_file="/home/vipink/Documents/FANTOM6/data/annotation/GRCh38-ELS.bed"

DNA_read_folder="/home/vipink/Documents/FANTOM6/data/RADICL_data/IPSC/replicate1/raw/DNA/chr/"
#%%
chr_files = os.listdir(f"{DNA_read_folder}")
trx_df = pd.read_csv(f"{annotation_file}",header=None,delimiter="\t")
trx_df.columns = ['chrom','start','end',
                  'ID','score','strand',
                  'start_b','end_b','sym',
                  'exon_count','exon_length','exon_start']
cluster_trx = bf.merge(trx_df)
enh_df = pd.read_csv(f"{enh_file}",header=None,delimiter="\t")
enh_df.columns = ['chrom','start','end','ID','ID2','type']

#%%
trx_dfs=[]
enh_dfs=[]
#%%
chr_set = list(map(lambda n: re.split('_|\\.',n)[3],os.listdir(f"{DNA_read_folder}")))
#%%
for chromo in chr_set:
    print(chromo)
    #%%
    read_df = pd.read_csv(f"{DNA_read_folder}RADICL_iPSC_DNA_{chromo}.bed",header=None,delimiter="\t")
    read_df.columns = ['chrom','start','end','ID','score']
    #%%
    trx_dfs.append(bf.count_overlaps(read_df,cluster_trx.query('chrom == @chromo'))
      .assign(genic=lambda df_:(df_.loc[:,'count']
                        .where(df_.loc[:,'count'].gt(0),'out')
                        .mask(df_.loc[:,'count'].gt(0),'in')))
      .groupby(['chrom','genic'])['count']
      .agg('count')
      .reset_index()
      .assign(set='transcript'))
    enh_dfs.append(bf.count_overlaps(read_df,enh_df.query('chrom == @chromo'))
      .assign(genic=lambda df_:(df_.loc[:,'count']
                        .where(df_.loc[:,'count'].gt(0),'out')
                        .mask(df_.loc[:,'count'].gt(0),'in')))
      .groupby(['chrom','genic'])['count']
      .agg('count')
      .reset_index()
      .assign(set='enhancers')

    )

# %%
pd.concat(enh_dfs).groupby('genic').agg(count=('count','sum'))
# %%
