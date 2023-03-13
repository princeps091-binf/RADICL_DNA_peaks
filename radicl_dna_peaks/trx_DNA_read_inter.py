#%%
import pandas as pd
#47843486 total reads in RADICL
#33373460 reads intersecting transcript
#5636408 reads intersecting enhancers
trx_read_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/DNA_read_transcript_intersect.bed"
enh_read_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/result/DNA_peak_analysis/DNA_read_enhancer_intersect.bed"

DNA_read_file="/home/vipink/Documents/FANTOM6/data/RADICL_data/raw/DNA/RADICL_iPSC_DNA.bed"
#%%
read_df = pd.read_csv(enh_read_file,header=None,delimiter="\t")

# %%
read_df.iloc[:,[3]].nunique()
# %%
