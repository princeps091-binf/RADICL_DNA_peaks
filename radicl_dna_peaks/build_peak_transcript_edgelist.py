#%%
import igraph as ig
import pandas as pd

#%%
peak_inter_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/peak_DNA_read_inter_tbl.tsv"
transcript_inter_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/transcript_to_peak_reads_map.tsv"

#%%
transcript_read_inter_df = pd.read_csv(transcript_inter_file,delimiter="\t")
peak_read_inter_tbl = pd.read_csv(peak_inter_file,delimiter="\t")

# %%
peak_transcript_map_tbl = (peak_read_inter_tbl
 .merge(transcript_read_inter_df,how='inner',left_on='read_ID',right_on='read_ID'))

#%%
peak_to_transcript_edge_list = (peak_transcript_map_tbl
 .loc[:,['chrom_x','start_x','end_x','ID','read_ID']]
 .assign(peak_ID=lambda df_:df_.chrom_x+ "_" + df_.start_x.astype(str)+ "_" + df_.end_x.astype(str))
 .loc[:,['peak_ID','ID','read_ID']]
 .groupby(by=['peak_ID','ID'])
 .count()
 .sort_values( by ='read_ID')
 .reset_index())
# %%
peak_trx_graph=ig.Graph.DataFrame(peak_to_transcript_edge_list,directed=False,use_vids=False)
# %%
# %%
peak_trx_graph.vs["type"] = pd.Series(peak_trx_graph.vs["name"]).isin(peak_to_transcript_edge_list.peak_ID.unique()).to_list()

# %%
peak_graph = peak_trx_graph.bipartite_projection(which=True)
# %%
peak_cmpnt = peak_graph.connected_components()
# %%
