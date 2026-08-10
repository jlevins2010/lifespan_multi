#!/usr/bin/env python
# coding: utf-8

# ## Further cleaning of the integrated object
# Cleaning and annotation of the multiVI object. It also subsets to just paired cells. Run on CUDA.

# ### Import Libraries

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata as ad
import scanpy as sc
import matplotlib as mpl
import os
import os.path
from IPython.display import display
import re


# In[2]:


colors = {"DCT_CNT": "black",
               "Endothelium": "#7ae031",
               "Podocyte": "#b56a07", 
               "Stroma": "#794b82",
               "NPC": "#ff8000", 
               "PT": "#ff00d4", 
               "Int": "#698cff",
               "IC": "#191985", 
               "PEC": "#ff0011", 
               "LOH": "#235e00",
               "Immune Cells": '#757575',
              "iPT":"#f78e52"
         }


# In[ ]:


sc.set_figure_params(
    dpi=80,  
    dpi_save=300,  
    frameon=False, 
    vector_friendly=True, 
    format="pdf",  
)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['pdf.fonttype'] = 42 
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.size'] = 7
plt.rcParams['axes.titlesize'] = 8
plt.rcParams['axes.labelsize'] = 7
plt.rcParams['xtick.labelsize'] = 6
plt.rcParams['ytick.labelsize'] = 6

sc.settings.figdir = "/home/levinsj/MultiOme/Analysis/pythonNotebooks/pythonPlots"


# ### scanpy variables

# In[3]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2
res = 2.0
n_features = 250000


# ### SCVI import

# In[4]:


batch_size = 32
max_epochs = 2000
earlyStopping = True
learning_rate = 0.0005
scvi_layers = 3
scvi_latent = 30
fraction_expressed = 0.01

multiVI_LATENT_KEY = "X_multivi"
multiVI_CLUSTERS_KEY = "clusters_multivi"


# ### SCANPY Settings

# In[5]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2


# In[6]:


def keep_matching_elements(my_list, adata_var_names):
  return [element for element in my_list if element in adata_var_names]


# In[7]:


demographics = pd.read_csv('/home/levinsj/MultiOme/Sample_demographics.csv') 
demographics.index = demographics["SampleName"]


# In[8]:


inputFile = "/home/levinsj/MultiOme/MergedObjects/merged_postMultiVI_V3_human_filtered_final.h5ad"
outputFile = "//home/levinsj/MultiOme/MergedObjects/merged_postMultiVI_annotated_final_human.h5ad"
output_barcodes_with_annotations = "/home/levinsj/MultiOme/Barcodes/final_annotated/all_annotated.csv"
output_barcodes_with_annotations_paired = "/home/levinsj/MultiOme/Barcodes/final_annotated/all_annotated_paired.csv"
anno_barcodes = "/home/levinsj/MultiOme/Barcodes/snapATAC_Barcodes/type_celltype_anno_barcodes.csv"


outputSCVIfile_genes_only = "/home/levinsj/MultiOme/MergedObjects/merged_postMultiVI_annotated_onlyGEX_human.h5ad"
outputSCVIfile_ATAC_only = "/home/levinsj/MultiOme/MergedObjects/merged_postMultiVI_annotated_onlyATAC_human.h5ad"
outputSCVIfile_paired_only = "/home/levinsj/MultiOme/MergedObjects/merged_postMultiVI_annotated_onlyPaired_human.h5ad"


# In[9]:


adata = sc.read_h5ad(inputFile)
print(adata)
sc.pp.neighbors(adata, use_rep=multiVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=min_distUMAP)
sc.tl.leiden(adata, key_added=multiVI_CLUSTERS_KEY, resolution=res)


# In[10]:


adata = adata[~adata.obs[multiVI_CLUSTERS_KEY].astype(str).isin(["52"])].copy()

cell_identities = {
    '0':'DCT_CNT','1':'LOH','2':'LOH','3':'NPC','4':'Endothelium','5':'PT','6':'PT',
    '7':'DCT_CNT','8':'LOH','9':'LOH','10':'Podocyte','11':'PT','12':'PT',
    '13':'PT','14':'Stroma', '15':'DCT_CNT', '16':'LOH', '17':'Stroma',
    '18':'Stroma','19':'Endothelium','20':'PT','21':'LOH', '22':'iPT',
    '23':'NPC','24':'Stroma','25':'PEC','26':'IC','27':'Stroma', '28':'DCT_CNT',
    '29':'Stroma', '30':'Immune Cells', '31':'LOH','32':'DCT_CNT', '33': 'Int', '34':'LOH',
    '35':'Immune Cells', '36':'LOH','37':'PT', '38':'Int','39':'LOH','40': 'LOH', '41':'Stroma',
    '42':'IC','43':'Endothelium','44':'DCT_CNT', '45':'Endothelium','46':'Int','47':'LOH','48':'DCT_CNT',
    '49': 'DCT_CNT', '50':'Endothelium','51':'Podocyte','53':'LOH','54':'DCT_CNT'
}

adata.obs["cellType"] = (
    adata.obs[multiVI_CLUSTERS_KEY]
    .astype(str)
    .map(cell_identities)
    .astype('category')
)


# In[11]:


ax = sc.pl.umap(
    adata,
    color="cellType",
    legend_loc='on data',
    legend_fontsize=7,
    legend_fontoutline=2,
    palette=colors,
    show=False                
)

fig = ax.figure

fig.savefig(
    '/home/levinsj/MultiOme/Analysis/pythonNotebooks/pythonPlots/UMAP_cellType_plot.pdf', 
    format='pdf', 
    bbox_inches='tight'
)

plt.close(fig)


# In[ ]:


ax = sc.pl.umap(
    adata,
    color="type",
    legend_loc='on data',
    legend_fontsize=7,
    legend_fontoutline=2,
    show=False                
)

fig = ax.figure

fig.savefig(
    '/home/levinsj/MultiOme/Analysis/pythonNotebooks/pythonPlots/UMAP_type_plot.pdf', 
    format='pdf', 
    bbox_inches='tight'
)

plt.close(fig)


# In[13]:


var_names = ['UNCX', 'ITGA8', 'LHX1', 'JAG1', 'PTPRO', 'MAFB',"CHF","VCAM1","HAVCR1",'CUBN', 'SLC13A1', 'SLC12A1', 'UMOD', 'SLC12A3', 'GATA3', 'SLC26A4', 'ATP6V1B1', 'PLVAP', 'EMCN', 'COL1A1', 'COL3A1', 'DOCK2', 'PTPRC']
order = ["NPC","Int","Podocyte","PEC", "iPT","PT","LOH","DCT_CNT","IC","Endothelium","Stroma","Immune Cells"]


# In[14]:


dotplot_obj = sc.pl.dotplot(
    adata,
    var_names=keep_matching_elements(var_names, adata.var_names),
    groupby='cellType',
    categories_order=order,
    standard_scale='var',
    cmap='Blues',
    dot_max=0.6,
    dot_min=0.0,
    figsize=(8, 6),
    dendrogram=False,
    log=True
)


# In[15]:


dotplot_obj = sc.pl.dotplot(
    adata,
    var_names=keep_matching_elements(var_names, adata.var_names),
    groupby='cellType',
    categories_order=order,
    standard_scale='var',
    cmap='Blues',
    dot_max=0.6,
    dot_min=0.0,
    figsize=(8, 6),
    dendrogram=False,
    log=True,
    show=False 
)

fig = plt.gcf()

fig.savefig(
    'dotPlot_allCells.pdf', 
    format='pdf', 
    bbox_inches='tight'
)

plt.close(fig)


# ## Demultiplex the adult samples

# In[16]:


import pandas as pd
import re

subdivision_configs = {
    'HK2801_HK3032': '/home/levinsj/MultiOme/Barcodes/DemultiPlex/HK2801_HK3032/final.csv',
    'HK2804_HK2826': '/home/levinsj/MultiOme/Barcodes/DemultiPlex/HK2804_HK2826/final.csv',
    'HK2823_HK3101': '/home/levinsj/MultiOme/Barcodes/DemultiPlex/HK2823_HK3101/final.csv',
}

extracted_names = adata.obs_names.str.extract(r'([ACGT]{16}(?:-\d)?)')[0]
print(extracted_names[0])
adata.obs_names = extracted_names.fillna(adata.obs_names.to_series()).astype(str)
print(adata.obs_names[0])


# In[17]:


adata.obs["sample_final"] = adata.obs["sample"].astype(str)

# 3. Process each CSV
for parent_sample_to_subdivide, csv_path in subdivision_configs.items():
    print(f"\nProcessing: {parent_sample_to_subdivide}")

    # Read CSV, set the barcode as the index
    sub_barcode_df = pd.read_csv(csv_path, index_col=0)
    sub_barcode_df.index = sub_barcode_df.index.astype(str)

    parent_mask = adata.obs["sample"] == parent_sample_to_subdivide
    cells_in_csv_mask = parent_mask & adata.obs.index.isin(sub_barcode_df.index)

    match_count = cells_in_csv_mask.sum()
    print(f"Found {match_count} matching cells in AnnData for this file.")

    if match_count > 0:
        # Pull matching barcodes and map them to the 'singlet' column
        matching_barcodes = adata.obs[cells_in_csv_mask].index
        adata.obs.loc[cells_in_csv_mask, 'sample_final'] = matching_barcodes.map(sub_barcode_df['singlet'])

print("\n--- Final Value Counts ---")
print(adata.obs["sample_final"].value_counts())


# ### update with a numeric gestational age for all samples

# In[18]:


demographics = pd.read_csv('/home/levinsj/MultiOme/Sample_demographics.csv') 
demographics.index = demographics["SampleName"]

keys = demographics['SampleName']
age = demographics['Age']
gAge = demographics['gAge']

# Use zip to pair them up and then convert to a dictionary
age_dict = dict(zip(keys, age))
gAge_dict = dict(zip(keys, gAge))

age_dict['HK2801'] = '45'
age_dict['HK3032'] = '61'
age_dict['HK2826'] = '68'
age_dict['HK2804'] = '49'
age_dict['HK3101'] = '48'
age_dict['HK2823'] = '20'
gAge_dict['HK2801'] = 'Not documented'
gAge_dict['HK3032'] = 'Not documented'
gAge_dict['HK2826'] = 'Not documented'
gAge_dict['HK2804'] = 'Not documented'
gAge_dict['HK3101'] = 'Not documented'
gAge_dict['HK2823'] = 'Not documented'

adata.X = adata.layers["counts"]
adata.obs["age"] = adata.obs["sample_final"]
adata.obs["gestational_age"] = adata.obs["sample_final"]

adata.obs["age"] = adata.obs["age"].map(age_dict).astype('category')
values_to_replace = ['45/61', '49/68', '20/48']
adata.obs["age"] = adata.obs["age"].replace(values_to_replace, np.nan)

non_nan_age_mask = adata.obs['age'].notna()
adata_cleaned = adata[non_nan_age_mask, :].copy()
adata_cleaned.obs['age'] = adata_cleaned.obs['age'].astype(float)

adata_cleaned.obs["gestational_age"] = adata_cleaned.obs["gestational_age"].map(gAge_dict).astype('category')
adata_cleaned.obs["gestational_age"] = adata_cleaned.obs["gestational_age"].str.slice(0, 2)
adata_cleaned.obs["gestational_age"] = adata_cleaned.obs["gestational_age"].replace("No", np.nan)
adata_cleaned.obs['gestational_age'] = adata_cleaned.obs['gestational_age'].astype(float)
adata_cleaned.obs['gestational_age'] = 0-((39-adata_cleaned.obs['gestational_age'])/52)

has_gestational_age = adata_cleaned.obs['gestational_age'].notna()

adata_cleaned.obs.loc[has_gestational_age, 'age'] = adata_cleaned.obs.loc[has_gestational_age, 'gestational_age']

print(adata_cleaned.obs["sample_final"].value_counts())
print(adata_cleaned.obs["type"].value_counts())

adata = adata_cleaned


# ### save each modality separately

# In[19]:


adata_genes = adata[:, adata.var["modality"] == "0"]
print(adata_genes)

adata_ATAC = adata[:, adata.var["modality"] == "1"]
print(adata_ATAC)

adata_paired = adata[adata.obs["modality"] == "paired"]
print(adata_paired)


# ### save files

# In[20]:


adata.write(outputFile)
adata_ATAC.write(outputSCVIfile_ATAC_only)
adata_genes.write(outputSCVIfile_genes_only)
adata_paired.write(outputSCVIfile_paired_only)


# ### update demographics file

# In[21]:


shared_col = demographics.columns.get_loc('sharedBarcodes_finalCells')

for index, item in enumerate(demographics['sampleRawH5Matricies_snRNA']):
    sampleName = demographics['SampleName'][index]
    if demographics['SampleName'][index] in adata_paired.obs["sample"].values:
        adata_subset = adata_paired[adata_paired.obs["sample"] == sampleName]
        demographics.iloc[index, shared_col] = adata_subset.n_obs
    else:
        demographics.iloc[index, shared_col] = 0

demographics.to_csv('/home/levinsj/MultiOme/Sample_demographics.csv', index=False)


# In[22]:


barcodes = adata.obs_names.to_list()
cell_types = adata.obs['cellType'].to_list()
sample = adata.obs["sample"].to_list()

data = {'barcode_orig': barcodes, 'cell_type': cell_types, 'sample': sample}

df = pd.DataFrame(data)
df['barcode'] = df['barcode_orig'].str.extract(r'[:_](.*)')

df.to_csv(output_barcodes_with_annotations, index=False)

barcodes = adata_paired.obs_names.to_list()
cell_types = adata_paired.obs['cellType'].to_list()
sample = adata_paired.obs["sample"].to_list()
data = {'barcode_orig': barcodes, 'cell_type': cell_types, 'sample': sample}

df_paired = pd.DataFrame(data)
df_paired['barcode'] = df_paired['barcode_orig'].str.extract(r'[:_](.*)')

df_paired.to_csv(output_barcodes_with_annotations_paired, index=False)


# In[23]:


adata_paired.obs['annoCellType'] = adata_paired.obs['type'].astype(str) + '_' + adata_paired.obs['cellType'].astype(str)
print(adata_paired.obs['annoCellType'].value_counts())
cell_types = adata_paired.obs['annoCellType'].to_list()

data = {'barcode_orig': barcodes, 'annoCellType': cell_types, 'sample': sample}

df_paired = pd.DataFrame(data)
df_paired['barcode'] = df_paired['barcode_orig'].str.extract(r'[:_](.*)')

df_paired.to_csv(anno_barcodes, index=False)


# ### calculate DEGs

# adata2 = adata_genes.copy()
# adata2.X = adata2.layers["counts"]
# sc.pp.normalize_total(adata2, target_sum=1e4)
# sc.pp.log1p(adata2)
# sc.tl.rank_genes_groups(adata2, multiVI_CLUSTERS_KEY, method='wilcoxon')
# sc.pl.rank_genes_groups(adata2, n_genes=25, sharey=False)
# leidenCheck = np.unique(adata2.obs["cellType"])
# sc.tl.rank_genes_groups(adata2, "cellType", method='wilcoxon')
# 
# for i in leidenCheck:
#     df = sc.get.rank_genes_groups_df(adata2, group = i)
#     fileName = "/home/levinsj/MultiOme/Analysis/DEGs/human/cellType_DEGs/"+str(i)+".csv"
#     df.to_csv(fileName)

# In[ ]:


ax = sc.pl.umap(
    adata,
    color="Age",
    legend_loc='on data',
    frameon=False,
    legend_fontsize=7,
    legend_fontoutline=2,
    show=False  
)

fig = ax.figure
fig.savefig('/home/levinsj/MultiOme/Analysis/pythonNotebooks/pythonPlots/UMAP_cellType_plot_age.pdf', format='pdf', bbox_inches='tight')


# In[ ]:




