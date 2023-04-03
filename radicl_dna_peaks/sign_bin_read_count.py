#%%
import pandas as pd
import os
import bioframe as bf
import re
#iPSC
# total readcount: 46041920
# in bin readcount: 43120168 -> 94%
# Neuron
# total readcount: 74004761
# in bin readcount: 70212907 -> 95%


#%%
chicane_anno_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/ChICANE_analysis_set18/interaction_ID.annotation.tsv"
chicane_interaction_sign_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/ChICANE_analysis_set18/CHICANE_significant_interaction_3cells.tsv"
DNA_read_folder="/home/vipink/Documents/FANTOM6/data/RADICL_data/Neuron/replicate1/raw/DNA/chr/"
#%%
chr_files = os.listdir(f"{DNA_read_folder}")
chicane_anno_df = pd.read_csv(chicane_anno_file,sep="\t",header=0)
chicane_inter_df = pd.read_csv(chicane_interaction_sign_file,sep="\t",header=0)
# %%
bin_to_ID_df = (chicane_anno_df
 .DNA_bin.str.split(pat="_",expand=True)
 .rename(columns={0:'chrom',1:"start",2:"end"})
 .assign(start=lambda df_:df_.start.astype('int64'),
         end=lambda df_:df_.end.astype('int64'))
 .assign(interaction_ID=chicane_anno_df.loc[:,'interaction_ID']))
bin_to_inter_sing_df = (chicane_inter_df
 .merge(bin_to_ID_df))
dna_bin_inter_df = (bin_to_inter_sing_df
 .query('sig_CHICANE_Neuron == "yes"')
 .groupby(['chrom','start','end'])
 .agg(inter_count=('interaction_ID','count'))
 .reset_index()
 .sort_values("inter_count"))
#%%
bin_dfs=[]
chr_set = list(map(lambda n: re.split('_|\\.',n)[2],os.listdir(f"{DNA_read_folder}")))
#%%
for chromo in chr_set:
    print(chromo)
    read_df = pd.read_csv(f"{DNA_read_folder}RADICL_DNA_{chromo}.bed",header=None,delimiter="\t")
    read_df.columns = ['chrom','start','end','ID','score']
    bin_dfs.append(bf.count_overlaps(read_df,(bin_to_inter_sing_df.query('chrom == @chromo')))
      .assign(genic=lambda df_:(df_.loc[:,'count']
                        .where(df_.loc[:,'count'].gt(0),'out')
                        .mask(df_.loc[:,'count'].gt(0),'in')))
      .groupby(['chrom','genic'])['count']
      .agg('count')
      .reset_index())

# %%
pd.concat(bin_dfs).groupby('genic').agg(read_count=('count','sum'))
# %%
