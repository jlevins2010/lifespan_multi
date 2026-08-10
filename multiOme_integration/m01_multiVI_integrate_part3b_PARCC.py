#!/usr/bin/env python
# coding: utf-8

# ## Further cleaning of the integrated object
# Cleaning and annotation of the multiVI object. It also subsets to just paired cells.

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


# ### scanpy variables

# In[2]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2
res = 2.0
n_features = 250000


# ### SCVI import

# In[3]:


batch_size = 64
max_epochs = 2000
earlyStopping = True
learning_rate = 0.0005
scvi_layers = 3
scvi_latent = 30
fraction_expressed = 0.01

multiVI_LATENT_KEY = "X_multivi"
multiVI_CLUSTERS_KEY = "clusters_multivi"


# ### SCANPY Settings

# In[4]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2


# ### demographics

# In[5]:


demographics = pd.read_csv('/vast/projects/chenyuli/nephrobase/levinsj/input_csv/Sample_demographics.csv') 
demographics.index = demographics["SampleName"]


# In[6]:


preSCVIfile = "/vast/projects/chenyuli/nephrobase/levinsj/adata/merged_preMultiVI_human.h5ad"
preSCVIfile_filtered = "/vast/projects/chenyuli/nephrobase/levinsj/adata/merged_preMultiVI_human_filtered.h5ad"
inputfile = "/vast/projects/chenyuli/nephrobase/levinsj/adata/merged_postMultiVI_V3_human_filtered.h5ad"
outputfile = "/vast/projects/chenyuli/nephrobase/levinsj/adata/merged_postMultiVI_V3_human_filtered_final.h5ad"


# In[7]:


adata = sc.read_h5ad(filename = inputfile)


# ### pulling full count data back to integrated object

# In[8]:


adata_original = sc.read_h5ad(preSCVIfile)


# In[9]:


obs_names_main = adata.obs_names
obs_names_original = adata_original.obs_names

# Find the intersection of the observation names
shared_obs_names = obs_names_main.intersection(obs_names_original)

print(f"Number of shared observations: {len(shared_obs_names)}")
print("-" * 30)

adata_original = adata_original[list(shared_obs_names), :].copy()

for key, matrix in adata.obsm.items():
    # Create a DataFrame from the matrix for easier indexing by obs_names
    temp_df = pd.DataFrame(matrix, index=adata.obs_names)

    # Select only the rows corresponding to shared_obs_names, ensuring correct order
    # .reindex will fill with NaNs for any missing cells (though they should all be present here)
    # and ensures the order matches adata_subset.obs_names
    aligned_matrix = temp_df.loc[adata.obs_names].values

    # Assign the aligned matrix to the .obsm of the subset object
    adata_original.obsm[key] = aligned_matrix
    print(f"Copied '{key}' from adata_reference to adata_subset.")

print("-" * 30)
print("adata_original .obsm keys (after transfer):", adata_original.obsm.keys())
print(adata_original)


# In[10]:


NPC_markers =  ["UNCX","SIX2","SIX1","CITED1"]
PT_markers  =  ["CUBN","APOE"]
Int_markers =  ["JAG1","LHX1"]
PEC_markers =  ["CFH","WT1","PTPRO"]
DCT_markers =  ["WNK1","AQP3","SLC12A3","TFAP2A"]
LOH_markers =  ["UMOD","AQP1","SLC12A1"]
Podo_markers = ["MAFB","NPHS2","NPHS1"]
Stro_markers = ["COL3A1","COL1A1",'PDGFRA']
Endo_markers = ["EFGL7", "PLVAP","EMCN"]
CNT_markers = ["TRPV5","CALB1","VDR","KL"]
IC_markers = ["SLC4A1",'ADGRF5','SLC26A7','SLC26A4','ATP6V1B1']
UB_markers = ["RET",'GATA3']


def keep_matching_elements(my_list, adata_var_names):
  return [element for element in my_list if element in adata_var_names]

sc.pl.umap(adata_original, color = keep_matching_elements(NPC_markers, adata_original.var_names), color_map = "viridis_r")
sc.pl.umap(adata_original, color = keep_matching_elements(Int_markers, adata_original.var_names), color_map = "viridis_r")
sc.pl.umap(adata_original, color = keep_matching_elements(Stro_markers, adata_original.var_names), color_map = "viridis_r")
sc.pl.umap(adata_original, color = keep_matching_elements(Endo_markers, adata_original.var_names), color_map = "viridis_r")
sc.pl.umap(adata_original, color = keep_matching_elements(PT_markers, adata_original.var_names), color_map = "viridis_r")
sc.pl.umap(adata_original, color = keep_matching_elements(IC_markers, adata_original.var_names), color_map = "viridis_r")


# In[11]:


print(adata_original.obs["sample"].value_counts())


# In[12]:


adata_original.write_h5ad(filename = outputfile)

