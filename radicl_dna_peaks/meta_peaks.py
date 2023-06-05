#%%
import pandas as pd
import numpy as np
import hdbscan
import networkx as nx
import bioframe as bf
import hvplot.pandas
import holoviews as hv

import multiprocessing
#%%
peak_file = "/home/vipink/Documents/FANTOM6/RADICL_hotspot/workflow/data/results/DNA/IPSC_replicate1/IPSC_replicate1_all_peak.bed"
read_stem = "/home/vipink/Documents/FANTOM6/RADICL_hotspot/workflow/data/processed/DNA/chr/IPSC_replicate1_clean_DNA_"
bdg_stem = "/home/vipink/Documents/FANTOM6/RADICL_hotspot/workflow/data/results/DNA/IPSC_replicate1/bdg/"
#%%
peak_df = pd.read_csv(peak_file,header=None,delimiter="\t")
peak_df.columns = ['chrom','start','end','ID','score',"strand","enrichment","pvalue","qvalue","peak"]
# %%
dist_df = bf.closest(peak_df.drop_duplicates())
# %%
plot = dist_df.assign(ld = lambda df_:list(np.log10(df_.distance))).hvplot.kde('ld')
hvplot.save(plot, 'test.html')



#%%
chromo = "chr1"
chr_peak_df = (peak_df
               .query('chrom == @chromo')
               .assign(width = lambda df_:df_.end-df_.start)
               .assign(mid = lambda df_:df_.start + df_.width/2))

#%%
# HDBSCAN on peak mid-points

clusterer = hdbscan.HDBSCAN(min_cluster_size=2,
                            cluster_selection_epsilon=25,
                            metric='euclidean')
clusterer.fit(chr_peak_df.loc[:,['mid']])
#%%
egde_list = clusterer.condensed_tree_.to_pandas()
#%%
# identify likely peak-aggregates (meta-peaks -> maximum entropy)

def cluster_nr(i):
    tmp_cl = (clusterer.single_linkage_tree_
                        .get_clusters(1/i,min_cluster_size=2))
    return sum(np.unique(tmp_cl) >= 0)
#%%
lambda_val = np.unique(egde_list.lambda_val.sort_values(ascending=True))[::-1][1:-1]

with multiprocessing.Pool(processes=5) as pool:
        # Using map_async method to perform square operation on all numbers parallely
        result_cln = pool.map(cluster_nr ,lambda_val)        
lambda_cl_df = pd.DataFrame({'distance':(1/np.array(lambda_val)),
                             'cl':result_cln})
#%%
lambda_cl_df.assign(ld=lambda df_:np.log10(df_.distance)).hvplot.line(x='ld',y='cl')
lambda_peak = 1/lambda_cl_df.sort_values(['cl','distance'],ascending=[False,True]).distance.iloc[0]
#%%
chr_cl_peak_df = (chr_peak_df
                    .assign(cl=clusterer.single_linkage_tree_
                                            .get_clusters(1/lambda_peak,min_cluster_size=2))
                    .groupby('cl')
                    .agg( chrom = ('chrom','first'),
                          start = ('start','min'),
                          end = ('end','max'),
                          npeak = ('start','count'))
                    .assign(width=lambda df_:df_.end-df_.start)
                    .reset_index()
                    .sort_values('npeak')
                    .query('cl > -1'))
chr_isolated_peaks_df = (chr_peak_df
                        .assign(cl=clusterer.single_linkage_tree_
                                            .get_clusters(1/lambda_peak,min_cluster_size=2))   
                        .query('cl < 0')
                        .loc[:,['cl','chrom','start','end','width']]
                        .assign(npeak=1))
chr_processed_peak_df = pd.concat([chr_cl_peak_df,chr_isolated_peaks_df])
#%%
# subset corresponding raw RADICL-reads to then re-apply hdbscan to unpack peak architecture 
read_file = f"{read_stem}{chromo}.bed"
DNA_read_df = (pd.read_csv(read_file,delimiter = "\t",header=None)
               .rename(columns={
                     0:'chrom',
                     1:'start',
                     2:'end',
                     3:'ID',
                     4:'score',
                     5:'strand'
               }))
#%%
chr_bdg_df = (pd.read_csv(f"{bdg_stem}{chromo}.bdg",delimiter = "\t",header=None)
              .rename(columns={
                    0:'chrom',
                    1:'start',
                    2:'end',
                    3:'qpois'
              }))
#%%
(bf.count_overlaps(chr_processed_peak_df,DNA_read_df)
 .assign(rate=lambda df_:df_.loc[:,'count']/df_.width)
 .sort_values('count')
 .loc[:,'count'].sum())

# will need to develop way to unpack the hierarchical configuration of high-density regions (read-level hdbscan + bedgraph)
# %%
(bf.count_overlaps(chr_peak_df
 .loc[:,['chrom','start','end','width']],DNA_read_df)
 .assign(rate=lambda df_:df_.loc[:,'count']/df_.width)
 .sort_values('rate')
 .loc[:,'count'].sum())
# %%
qpois_track_df = bf.overlap(chr_bdg_df,chr_processed_peak_df).query("chrom_.notnull()")
# %%
plot = qpois_track_df.hvplot.errorbars(xerr1='start',xerr2='end',y='qpois')
# %%
plot = hv.Segments(qpois_track_df.query('start > 228608516').query('end < 228646431'),["start",'qpois','end','qpois'])
plot = plot.opts(line_width=50)

hvplot.save(plot, 'test.html')


# %%
