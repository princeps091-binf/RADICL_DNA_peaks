#%%
import pandas as pd
import bioframe as bf
import altair as alt
#%%
neuron_peak_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/Neuron/replicate1/results/MACS/peaks/RADICL_DNA_tot_peaks.bed"
ipsc_peak_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/IPSC/replicate1/results_beta/MACS/peaks/RADICL_DNA_tot_peaks.bed"
black_list_file = "/home/vipink/Documents/FANTOM6/data/annotation/hg38-blacklist.v2.bed"

#%%
neuron_peak_df=pd.read_csv(neuron_peak_file,sep="\t",header=None)
neuron_peak_df.columns = ['chrom','start','end','ID','scoreA',"sym",'scoreB','scoreC','scoreD','scoreE']
ipsc_peak_df=pd.read_csv(ipsc_peak_file,sep="\t",header=None)
ipsc_peak_df.columns = ['chrom','start','end','ID','scoreA',"sym",'scoreB','scoreC','scoreD','scoreE']
black_list_df=pd.read_csv(black_list_file,sep="\t",header=None)
black_list_df.columns = ['chrom','start','end','label']

# %%
bf.count_overlaps(neuron_peak_df,black_list_df).query('count > 0')
# %%
