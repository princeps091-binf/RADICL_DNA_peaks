#%% 
import pandas as pd
import bioframe as bf
import altair as alt
#%%
chicane_anno_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/ChICANE_analysis_set18/interaction_ID.annotation.tsv"
chicane_interaction_sign_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/ChICANE_analysis_set18/CHICANE_significant_interaction_3cells.tsv"
neuron_peak_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/Neuron/replicate1/results/MACS/peaks/RADICL_DNA_tot_peaks.bed"
ipsc_peak_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/IPSC/replicate1/results_beta/MACS/peaks/RADICL_DNA_tot_peaks.bed"
#%%
chicane_anno_df = pd.read_csv(chicane_anno_file,sep="\t",header=0)
chicane_inter_df = pd.read_csv(chicane_interaction_sign_file,sep="\t",header=0)
neuron_peak_df=pd.read_csv(neuron_peak_file,sep="\t",header=None)
neuron_peak_df.columns = ['chrom','start','end','ID','scoreA',"sym",'scoreB','scoreC','scoreD','scoreE']
ipsc_peak_df=pd.read_csv(ipsc_peak_file,sep="\t",header=None)
ipsc_peak_df.columns = ['chrom','start','end','ID','scoreA',"sym",'scoreB','scoreC','scoreD','scoreE']

# %%
bin_to_ID_df = (chicane_anno_df
 .DNA_bin.str.split(pat="_",expand=True)
 .rename(columns={0:'chrom',1:"start",2:"end"})
 .assign(start=lambda df_:df_.start.astype('int64'),
         end=lambda df_:df_.end.astype('int64'))
 .assign(interaction_ID=chicane_anno_df.loc[:,'interaction_ID']))
# %%
bin_to_inter_sing_df = (chicane_inter_df
 .merge(bin_to_ID_df))
# %%
dna_bin_inter_df = (bin_to_inter_sing_df
 .query('sig_CHICANE_iPSC == "yes"')
 .groupby(['chrom','start','end'])
 .agg(inter_count=('interaction_ID','count'))
 .reset_index()
 .sort_values("inter_count"))
# %%
dna_bin_inter_df = bf.count_overlaps(dna_bin_inter_df,ipsc_peak_df)
# %%
alt.data_transformers.disable_max_rows()

base = (alt.Chart(dna_bin_inter_df.query('count > 0'))
.mark_point(
    size=1,
    filled=True,
    opacity=1
)
.encode(
    alt.X("count:Q"),
    alt.Y('inter_count:Q')))

base + base.transform_regression('count', 'inter_count',method="poly",order=1).mark_line(size=4)
# %%
(bf.count_overlaps(ipsc_peak_df,dna_bin_inter_df)
 .query('count > 0')
 .sort_values('count'))
# %%
