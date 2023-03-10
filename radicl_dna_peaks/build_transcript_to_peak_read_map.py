#%%
import bioframe as bf
import pandas as pd
import subprocess
#%%
annotation_file = "/home/vipink/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed"
rna_radicl_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/RADICL_peak_RNA.bed"
tmp_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/tmp/tmp_res.txt"
out_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/transcript_to_peak_reads_map.tsv"
#%%
rna_df = pd.read_csv(rna_radicl_file,header=None,delimiter="\t")
rna_df.columns = ['chrom','start','end','read_ID','score','strand']

# %%
transcript_annotation_df = pd.read_csv(annotation_file,header=None,delimiter="\t")

transcript_annotation_df.columns = ['chrom','start','end',
                                    'ID','score','strand',
                                    'start_b','end_b','sym',
                                    'exon_count','exon_length','exon_start']

# %%
transcript_read_inter_idx = bf.overlap(transcript_annotation_df, rna_df, 
                                       on=['strand'],
                                       how='inner', 
                                       suffixes=('_1','_2'),
                                       return_index=True,
                                       return_input=False)
# %%
transcript_read_inter_tbl = (transcript_annotation_df.loc[transcript_read_inter_idx.index_1,["chrom","start","end","strand","ID"]].reset_index(drop=True)
 .assign(read_ID=rna_df.loc[transcript_read_inter_idx.index_2,["read_ID"]].reset_index(drop=True)))

# %%
transcript_read_inter_tbl.to_csv(out_file,
                           sep="\t",
                           header=True,
                           index=False)
