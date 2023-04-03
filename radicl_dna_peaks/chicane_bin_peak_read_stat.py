#%%
import pandas as pd
import bioframe as bf
import scipy.stats as stats
import altair as alt
#iPSC
# total readcount: 46041920
# in bin readcount: 43120168 -> 94%
# Neuron
# total readcount: 74004761
# in bin readcount: 70212907 -> 95%

#%%
peak_read_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/python_project_result/DNA_peak_analysis/Neuron_rep1_peak_DNA_read_inter_tbl.tsv"
chicane_anno_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/ChICANE_analysis_set18/interaction_ID.annotation.tsv"
chicane_interaction_sign_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/ChICANE_analysis_set18/CHICANE_significant_interaction_3cells.tsv"
#%%
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
peak_read_tbl = pd.read_csv(peak_read_file,delimiter="\t")
#%%
peak_read_in_trx = (bf.count_overlaps(peak_read_tbl,dna_bin_inter_df)
 .query('count>0')
 .read_ID.size)
peak_read_tot = peak_read_tbl.read_ID.size
# %%
obs_mat = [[peak_read_in_trx,70212907],[peak_read_tot,74004761]]
chi_result = stats.chi2_contingency(obs_mat)
# %%
(obs_mat - chi_result.expected_freq)[0,:]

# %%
tmp_dat=[peak_read_in_trx/peak_read_tot,
         (peak_read_tot-peak_read_in_trx)/peak_read_tot,
         70212907/74004761,
         (74004761-70212907)/74004761]
tmp_dat = (pd.DataFrame({
    'count':tmp_dat,
    'in_out':['in','out','in','out'],
    'set':['peak','peak','total','total']
}))
# %%
(alt.Chart(tmp_dat,
           title=alt.TitleParams("significant bins"))
 .mark_bar().encode(
    x="set",
    y="count",
    color="in_out"
))

# %%
