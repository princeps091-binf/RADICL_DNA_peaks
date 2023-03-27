#%%
import pandas as pd
import numpy as np
import altair as alt

#%%
peak_to_read_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/peak_DNA_read_inter_tbl.tsv"
transcript_to_read_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/transcript_to_peak_reads_map.tsv"
peak_to_transcript_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/peak_to_transcript_tbl.tsv"
#%%
peak_to_read_df=pd.read_csv(peak_to_read_file,sep='\t')
transcript_to_read_df=pd.read_csv(transcript_to_read_file,sep='\t')
# %%
peak_to_transcript_df=(peak_to_read_df
 .merge(transcript_to_read_df,
        how='left',on='read_ID',
        suffixes=('_peak','_transcript')))

#%%
peak_to_transcript_df = (peak_to_transcript_df
 .assign(genic=lambda df_:(df_.strand
                        .where(df_.strand.isnull(),'genic')
                        .mask(df_.strand.isnull(),'intergenic'))))

#%%
(peak_to_transcript_df
 .to_csv(peak_to_transcript_file,
         sep='\t',
         index=False))

