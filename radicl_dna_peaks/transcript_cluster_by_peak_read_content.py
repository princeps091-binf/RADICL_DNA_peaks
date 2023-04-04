#%%
import pandas as pd
import bioframe as bf

#%%
transcript_to_read_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/python_project_result/DNA_peak_analysis/transcript_to_peak_reads_map.tsv"

#%%
transcript_to_read_df=pd.read_csv(transcript_to_read_file,sep='\t')

# %%
read_to_transcript_set_df = (transcript_to_read_df
 .groupby("read_ID")
 .agg(transcript_set = ('ID',lambda x: ','.join(sorted(list(x)))))
 .reset_index()
 )
# %%
(read_to_transcript_set_df
 .groupby('transcript_set')
 .agg(read_count = ('read_ID','count'))
 .reset_index()
 .sort_values('read_count')
)
# %%
