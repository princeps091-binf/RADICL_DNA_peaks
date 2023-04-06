#%%
import ibis
from ibis import _
#%%
RADICL_read_file = "/home/vipink/Documents/FANTOM6/data/RADICL_data/Neuron/replicate1/raw/15.GGCTAC.2.bed"

#%%
t = ibis.read_csv(RADICL_read_file,delim='\t')
# %%
t = t.relabel(dict(
    column00='RNA_chrom',
    column01='RNA_start',
    column02='RNA_end',
    column03='RNA_ID',
    column04='RNA_score',
    column05='RNA_strand',
    column06='DNA_chrom',
    column07='DNA_start',
    column08='DNA_end',
    column09='DNA_ID',
    column10='DNA_score',
    column11='DNA_strand',
    column12='RNA_label',
    column13='DNA_label',
    column14='inter_label'))
t.head().execute()
# %%
t.RNA_chrom.value_counts().order_by(_['RNA_chrom_count'].desc()).execute()
#%%
t_intra = t.filter(t.RNA_chrom == t.DNA_chrom).execute()
# %%
t_intra.inter_label.value_counts()
# %%
