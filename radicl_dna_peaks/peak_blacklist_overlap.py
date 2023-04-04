#%%
import pandas as pd
import bioframe as bf
import altair as alt
#%%
neuron_peak_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/Neuron/replicate1/results/MACS/peaks/RADICL_DNA_tot_peaks.bed"
ipsc_peak_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/IPSC/replicate1/results_beta/MACS/peaks/RADICL_DNA_tot_peaks.bed"
neuron_peak_read_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/python_project_result/DNA_peak_analysis/Neuron_rep1_peak_DNA_read_inter_tbl.tsv"
ipsc_peak_read_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/python_project_result/DNA_peak_analysis/IPSC_rep1_peak_DNA_read_inter_tbl.tsv"

black_list_file = "/home/vipink/Documents/FANTOM6/data/annotation/hg38-blacklist.v2.bed"

#%%
neuron_peak_df=pd.read_csv(neuron_peak_file,sep="\t",header=None)
neuron_peak_df.columns = ['chrom','start','end','ID','scoreA',"sym",'scoreB','scoreC','scoreD','scoreE']
ipsc_peak_df=pd.read_csv(ipsc_peak_file,sep="\t",header=None)
ipsc_peak_df.columns = ['chrom','start','end','ID','scoreA',"sym",'scoreB','scoreC','scoreD','scoreE']
#%%
black_list_df=pd.read_csv(black_list_file,sep="\t",header=None)
black_list_df.columns = ['chrom','start','end','label']
#%%
neuron_peak_read_df=pd.read_csv(neuron_peak_read_file,sep="\t",header=0)
ipsc_peak_read_df=pd.read_csv(ipsc_peak_read_file,sep="\t",header=0)

# %%
bf.count_overlaps(neuron_peak_df,black_list_df).query('count > 0')
#reads
#Neuron:918229/8400122 -> 10%
#iPSC:923796/1356429 -> 68 %
#peaks
#Neuron:9712/486995 -> 2%
#iPSC:4093/30400 -> 13 %

# %%
ipsc_peak_read_in_trx = (bf.count_overlaps(ipsc_peak_read_df,black_list_df)
 .query('count>0')
 .read_ID.size)
ipsc_peak_read_tot = ipsc_peak_read_df.read_ID.size

neuron_peak_read_in_trx = (bf.count_overlaps(neuron_peak_read_df,black_list_df)
 .query('count>0')
 .read_ID.size)
neuron_peak_read_tot = neuron_peak_read_df.read_ID.size


tmp_dat=[ipsc_peak_read_in_trx/ipsc_peak_read_tot,
         (ipsc_peak_read_tot-ipsc_peak_read_in_trx)/ipsc_peak_read_tot,
         neuron_peak_read_in_trx/neuron_peak_read_tot,
         (neuron_peak_read_tot-neuron_peak_read_in_trx)/neuron_peak_read_tot]
tmp_dat = (pd.DataFrame({
    'count':tmp_dat,
    'in_out':['in','out','in','out'],
    'set':['ipsc','ipsc','neuron','neuron']
}))
# %%
(alt.Chart(tmp_dat,
           title=alt.TitleParams("Black list regions"))
 .mark_bar().encode(
    x="set",
    y="count",
    color="in_out"
))

# %%
