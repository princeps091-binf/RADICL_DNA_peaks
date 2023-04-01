#%%
import pandas as pd
import bioframe as bf
import numpy as np
import altair as alt
#%%
rep1_peak_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/IPSC/replicate1/results_beta/MACS/peaks/RADICL_DNA_tot_peaks.bed"
rep2_peak_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/Neuron/replicate1/results/MACS/peaks/RADICL_DNA_tot_peaks.bed"

rep1_bdg_folder = "/home/vipink/Documents/FANTOM6/data/RADICL_data/IPSC/replicate1/results_beta/MACS/qvalue_track/"
rep2_bdg_folder = "/home/vipink/Documents/FANTOM6/data/RADICL_data/Neuron/replicate1/results/MACS/qvalue_track/"

# %%
rep1_peak_df = pd.read_csv(rep1_peak_file,header=None,delimiter="\t")
rep1_peak_df.columns = ['chrom','start','end','ID','score',"strand","enrichment","pvalue","qvalue","peak"]

rep2_peak_df = pd.read_csv(rep2_peak_file,header=None,delimiter="\t")
rep2_peak_df.columns = ['chrom','start','end','ID','score',"strand","enrichment","pvalue","qvalue","peak"]

#%%
# loop through each chromosome and assign qValue score for each peak
chr_set = rep2_peak_df.chrom.drop_duplicates().tolist()
#%%
dfs = []

for chromo in chr_set:
    print(chromo)
    chrom_rep1_peak = rep1_peak_df.query('chrom == @chromo')
    chrom_rep2_peak = rep2_peak_df.query('chrom == @chromo')

    rep1_bdg = pd.read_csv(f"{rep1_bdg_folder}RADICL_DNA_qval_{chromo}.bdg",header=None,delimiter="\t")
    rep1_bdg.columns = ['chrom','start','end','bdg']
    rep2_bdg = pd.read_csv(f"{rep2_bdg_folder}RADICL_DNA_qval_{chromo}.bdg",header=None,delimiter="\t")
    rep2_bdg.columns = ['chrom','start','end','bdg']
    score1 = (bf.overlap(pd.concat([chrom_rep1_peak,chrom_rep2_peak]),rep1_bdg)
     .groupby(['chrom','start','end'])
     .agg(qvalue1=('bdg_','max'))
     .reset_index()
     .assign(rank1=lambda df_:df_.qvalue1.rank(ascending=False)))
    score2 = (bf.overlap(pd.concat([chrom_rep1_peak,chrom_rep2_peak]),rep2_bdg)
     .groupby(['chrom','start','end'])
     .agg(qvalue2=('bdg_','max'))
     .reset_index()
     .assign(rank2=lambda df_:df_.qvalue2.rank(ascending=False)))
    dfs.append(score1.merge(score2))
# %%
alt.data_transformers.disable_max_rows()

scatter_chart = (alt.Chart(pd.concat(dfs).assign(tot_rank1=lambda df_:df_.qvalue1.rank(ascending=False),
                                 tot_rank2=lambda df_:df_.qvalue2.rank(ascending=False)))
.mark_point(size=1,opacity=0.5)
.encode(
    alt.X("tot_rank1:Q",scale=alt.Scale(type="log")),
    alt.Y('tot_rank2:Q',scale=alt.Scale(type="log")))
    .configure_point(
    filled=True
))

# %%

scatter_chart.save(f"/home/vipink/Documents/FANTOM6/litt/poster/img/neuron_vs_ipsc_rank.html", scale_factor=1.0)

# %%
