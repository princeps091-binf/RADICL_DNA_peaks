#%%
import pandas as pd
import bioframe as bf
import scipy.stats as stats
#%%
peak_read_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/peak_DNA_read_inter_tbl.tsv"
transcript_annotation_file="/home/vipink/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed"
transcript_cluster_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/transcript_cluster.bed"
cre_file="/home/vipink/Documents/FANTOM6/data/annotation/GRCh38-cCREs.bed"
enh_file="/home/vipink/Documents/FANTOM6/data/annotation/GRCh38-ELS.bed"
peak_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA/MACS/peaks/RADICL_DNA_tot_peaks.bed"
#%%
peak_read_tbl = pd.read_csv(peak_read_file,delimiter="\t")
transcript_annotation_df = pd.read_csv(transcript_annotation_file,header=None,delimiter="\t")
cre_df = pd.read_csv(cre_file,header=None,delimiter="\t")
peak_df = pd.read_csv(peak_file,header=None,delimiter="\t")
enh_df = pd.read_csv(enh_file,header=None,delimiter="\t")
transcript_cluster_df= pd.read_csv(transcript_cluster_file,header=0,delimiter="\t")
# %%
transcript_annotation_df.columns = ['chrom','start','end',
                                    'ID','score','strand',
                                    'start_b','end_b','sym',
                                    'exon_count','exon_length','exon_start']

enh_df.columns = ['chrom','start','end',
                  'ID1','ID2','label']

cre_df.columns = ['chrom','start','end',
                  'ID1','ID2','label']

peak_df.columns = ['chrom','start','end',
                  'ID','score','label',
                  'score_a','score_b','score_c','score_d']

# %%
peak_read_in_trx = (bf.count_overlaps(peak_read_tbl,enh_df)
 .query('count>0')
 .read_ID.size)
peak_read_tot = peak_read_tbl.read_ID.size
# %%
obs_mat = [[peak_read_in_trx,5636408],[peak_read_tot,47843486]]
chi_result = stats.chi2_contingency(obs_mat)
# %%
(obs_mat - chi_result.expected_freq)[0,:]
# %%
