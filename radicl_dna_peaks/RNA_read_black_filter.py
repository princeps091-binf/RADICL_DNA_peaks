#%%
import pandas as pd
import ibis 
import bioframe as bf
import os
#%%
black_list_file = "/home/vipink/Documents/FANTOM6/data/annotation/hg38-blacklist.v2.bed"
RNA_read_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/IPSC/replicate1/raw/RNA/RADICL_iPSC_RNA.bed"
transcript_annotation_file="/home/vipink/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed"
clean_rna_out="/home/vipink/Documents/FANTOM6/data/RADICL_data/IPSC/replicate1/processed/"

#%%
os.makedirs(f"{clean_rna_out}RNA/clean",exist_ok=True)

#%%
transcript_annotation_df = pd.read_csv(transcript_annotation_file,header=None,delimiter="\t")
transcript_annotation_df.columns = ['chrom','start','end',
                                    'ID','score','strand',
                                    'start_b','end_b','sym',
                                    'exon_count','exon_length','exon_start']

chr_set = transcript_annotation_df.chrom.drop_duplicates().to_numpy()
#%%
black_list_df=pd.read_csv(black_list_file,sep="\t",header=None)
black_list_df.columns = ['chrom','start','end','label']

#%%
read_df= ibis.read_csv(RNA_read_file,delim='\t')

# %%
read_df = read_df.relabel(dict(
    column0='chrom',
    column1='start',
    column2='end',
    column3='ID',
    column4='score',
    column5='strand'))

# %%
for chromo in chr_set:
    print(chromo)
    chr_read_df = read_df.filter(read_df.chrom == chromo).execute()
    
    clean_chr_read_df = bf.subtract(chr_read_df,black_list_df)
    
    clean_chr_read_df.to_csv(f"{clean_rna_out}RNA/clean/RADICL_RNA_{chromo}.bed",
                             sep="\t",
                             header=True,
                             index=False)
    del clean_chr_read_df
    del chr_read_df
# %%
