#%%
import bioframe as bf
import pandas as pd
import re
#%%
peak_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/Neuron/replicate1/results/MACS/peaks/RADICL_DNA_tot_peaks.bed"
DNA_read_folder = "/home/vipink/Documents/FANTOM6/data/RADICL_data/Neuron/replicate1/raw/DNA/chr/"
peak_read_inter_tbl_out_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/python_project_result/DNA_peak_analysis/Neuron_rep1_peak_DNA_read_inter_tbl.tsv"

#%%
peak_df = pd.read_csv(peak_file,header=None,delimiter="\t")
peak_df.columns = ['chrom','start','end','ID','score',"char","score_a","score_b","score_c","score_d"]

#%%
reads_dfs=[]
chr_set = list(map(lambda n: re.split('_|\\.',n)[2],os.listdir(f"{DNA_read_folder}")))
#%%
for chromo in chr_set:
    print(chromo)
    read_df = pd.read_csv(f"{DNA_read_folder}RADICL_DNA_{chromo}.bed",header=None,delimiter="\t")
    read_df.columns = ['chrom','start','end','ID','score']
    peak_read_inter_idx = bf.overlap(peak_df, read_df, how='inner', suffixes=('_1','_2'),return_index=True,return_input=False)
    reads_dfs.append(peak_df.loc[peak_read_inter_idx.index_1,["chrom","start","end"]].reset_index(drop=True)
                            .assign(read_ID=read_df.loc[peak_read_inter_idx.index_2,["ID"]].reset_index(drop=True)))

# %%
pd.concat(reads_dfs).to_csv(peak_read_inter_tbl_out_file,
                           sep="\t",
                           header=True,
                           index=False)
# %%
