#!/usr/bin/env python
# coding: utf-8

# ### Import Libraries

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata as ad
import scanpy as sc
import matplotlib as mpl
import torch
import os
import os.path
from IPython.display import display
import scvi
import muon as mu


# ### SCVI import

# In[2]:


scvi.__version__

scvi.settings.seed = 0
batch_size = 64
max_epochs = 2000
earlyStopping = True
learning_rate = 0.0001
scvi_layers = 3
scvi_latent = 30
fraction_expressed = 0.0005

torch.set_float32_matmul_precision('high')

multiVI_LATENT_KEY = "X_multivi"
multiVI_CLUSTERS_KEY = "clusters_multivi"


# ### SCANPY Settings

# In[3]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2


# In[ ]:


demographics = pd.read_csv('/vast/projects/chenyuli/nephrobase/levinsj/input_csv/Sample_demographics.csv') 
demographics.index = demographics["SampleName"]
display(demographics)


# ### Get shared barcodes

# In[5]:


snRNAinput = "/vast/projects/chenyuli/nephrobase/levinsj/adata/postSCVI_snRNA_annotated_human.h5ad"
snATACinput = "//vast/projects/chenyuli/nephrobase/levinsj/adata/postSCVI_ATAC_peaks_called.h5ad"

preSCVIfile = "/vast/projects/chenyuli/nephrobase/levinsj/adata/merged_preMultiVI_human.h5ad"
multiVImodel = "/vast/projects/chenyuli/nephrobase/levinsj/scvi_models/multiVI/"
postSCVIfile = "/vast/projects/chenyuli/nephrobase/levinsj/adata/merged_postMultiVI_human.h5ad"


# In[6]:


snRNA_adata = sc.read_h5ad(snRNAinput)
print(snRNA_adata.obs_names.duplicated().any())
snRNA_adata.var.modality = "Gene Expression"


# In[7]:


snATAC_adata = sc.read_h5ad(snATACinput)
print(snATAC_adata.obs_names.duplicated().any())
snATAC_adata.var.modality = "Peaks"
dictionary = pd.Series(demographics['PrepType'].values, index=demographics['SampleName']).to_dict()
snATAC_adata.obs["PrepType"] = snATAC_adata.obs['sample'].map(dictionary).astype('category')
snATAC_adata.obs["PrepType"] = snATAC_adata.obs["PrepType"].str.replace('snRNA/snATAC', 'snATAC')
snATAC_adata.obs["PrepType"] = snATAC_adata.obs["PrepType"].str.replace('MultiOme', 'multiOme')

dictionary = pd.Series(demographics['Site'].values, index=demographics['SampleName']).to_dict()
snATAC_adata.obs["study"] = snATAC_adata.obs['sample'].map(dictionary).astype('category')

dictionary = pd.Series(demographics['Type'].values, index=demographics['SampleName']).to_dict()
snATAC_adata.obs["type"] = snATAC_adata.obs['sample'].map(dictionary).astype('category')

print(snATAC_adata.obs["type"].value_counts())
print(snATAC_adata.obs["study"].value_counts())
print(snATAC_adata.obs["type"].value_counts())


# In[8]:


print(snRNA_adata)
print(snATAC_adata)


# In[9]:


sc.pl.highest_expr_genes(snRNA_adata)
sc.pl.highest_expr_genes(snATAC_adata)


# In[10]:


common_obs_names = snRNA_adata[snRNA_adata.obs["PrepType"] == "multiOme"].obs_names.intersection(snATAC_adata[snATAC_adata.obs["PrepType"] == "multiOme"].obs_names)

snRNA_paired = snRNA_adata[common_obs_names]
snATAC_paired = snATAC_adata[common_obs_names]

snRNA_unpaired = snRNA_adata[snRNA_adata.obs_names.difference(snRNA_paired.obs_names)]
snATAC_unpaired = snATAC_adata[snATAC_adata.obs_names.difference(snATAC_paired.obs_names)]


# In[11]:


adatas = {"snATAC": snATAC_paired,"snRNA": snRNA_paired}

adata_multiOme = ad.concat(
    [snRNA_paired, snATAC_paired],
    join="outer",  
    axis=1,        
    merge="first",  
    label="modality", 
)

print(adata_multiOme)

adata_multiOme.obs["class"] = "paired_multiOme"
snRNA_unpaired.obs["class"] = "snRNA_unpaired"
snATAC_unpaired.obs["class"] = "snATAC_unpaired"


# Combine the multiOme into one object

# In[12]:


res_anndata = adata_multiOme.copy()
modality_ann = ["paired"] * adata_multiOme.shape[0]
obs_names = list(adata_multiOme.obs.index.values)
modality_key = "modality"
obs_names_multiOme = obs_names


# In[13]:


adata_multiOme.obs.index = pd.Series(obs_names) + "_paired"
print(adata_multiOme.var.index.values)
print(adata_multiOme.obs.index.values)


# In[14]:


#remove .obsm so we can organize anndata (necessary for next step)
del adata_multiOme.obsm
del snRNA_unpaired.obsm
del snATAC_unpaired.obsm

adata_mvi = scvi.data.organize_multiome_anndatas(adata_multiOme, snRNA_unpaired, snATAC_unpaired)
sc.pl.highest_expr_genes(adata_mvi)

adata_mvi = adata_mvi[:, adata_mvi.var["modality"].argsort()].copy()
adata_mvi.obs['PrepType'] = adata_mvi.obs['PrepType'].fillna("multiOme")
adata_mvi.obs['total_counts'] = adata_mvi.obs['total_counts'].fillna(0)
adata_mvi.obs['pct_counts_mt'] = adata_mvi.obs['pct_counts_mt'].fillna(0)

sc.pl.highest_expr_genes(adata_mvi)


# Sort so that the gene expression modality is first

# In[15]:


print(adata_mvi.obs["PrepType"].value_counts())
print(adata_mvi.obs["study"].value_counts())
print(adata_mvi.obs["sample"].value_counts())
print(adata_mvi.var["modality"].value_counts())


# In[16]:


sc.pl.highest_expr_genes(adata_mvi)


# In[17]:


print(adata_mvi.var)


# In[18]:


adata = adata_mvi[:, adata_mvi.var["modality"].argsort()].copy()
print(adata.var)


# In[19]:


adata.layers["counts"] = adata.X.copy()

adata.write_h5ad(filename = preSCVIfile)


# In[20]:


sc.pl.highest_expr_genes(adata)

