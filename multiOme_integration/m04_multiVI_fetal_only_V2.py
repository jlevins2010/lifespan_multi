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


# In[2]:


colors = {"DCT_CNT": "black",
               "Endothelium": "#7ae031",
               "Podocyte": "#ad9c00", 
               "Stroma": "#794b82",
               "NPC": "#ff8000", 
               "PT": "#ff00d4", 
               "Int": "#698cff",
               "IC": "#191985", 
               "PEC": "#ff0011", 
               "LOH": "#235e00",
               "Immune Cells": '#757575',
         }


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


# In[5]:


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


# ### SCANPY Settings

# In[6]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2


# In[7]:


demographics = pd.read_csv('/home/levinsj/MultiOme/Sample_demographics.csv') 
demographics.index = demographics["SampleName"]


# In[8]:


inputFile_multiOme = "/home/levinsj/MultiOme/MergedObjects/merged_postMultiVI_annotated_final_human.h5ad"

fetal_only_multiOme = "/home/levinsj/MultiOme/MergedObjects/merged_postMultiVI_annotated_final_fetal.h5ad"
fetal_only_multiOme_nephron = "/home/levinsj/MultiOme/MergedObjects/merged_postMultiVI_annotated_final_fetal_nephron_only.h5ad"


# In[9]:


adata = sc.read_h5ad(inputFile_multiOme)
adata = adata[adata.obs["type"] == "Fetal"]

print(adata)

sc.pp.neighbors(adata, use_rep=multiVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=min_distUMAP)
adata.obs["cellType"][adata.obs["cellType"] == "iPT"] = "Int"
sc.pl.umap(adata, color = "cellType", legend_loc='on data', save="_cellType_fetal_only.pdf", show=True)
sc.pl.umap(adata, color = "RET",  cmap = "viridis_r", save="_RET_fetal_only.pdf", show=True)
sc.pl.umap(adata, color = "WNT9B", cmap = "viridis_r", save="_WNT9B_fetal_only.pdf", show=True)


adata.write_h5ad(filename = fetal_only_multiOme)

nephron_lineage = ["NPC","iPT","PT","LOH","DCT_CNT","PEC","Podocyte","Int"]
adata = adata[adata.obs["cellType"].isin(nephron_lineage)]

sc.pp.neighbors(adata, use_rep=multiVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=min_distUMAP)
sc.pl.umap(adata, color = "cellType", legend_loc='on data', save="_cellType_fetalNephron_only.pdf", show=True)
sc.pl.umap(adata, color = "RET",  cmap = "viridis_r", save="_RET_fetalNephron_only.pdf", show=True)
sc.pl.umap(adata, color = "WNT9B", cmap = "viridis_r", save="_WNT9B_fetalNephron_only.pdf", show=True)


sc.tl.leiden(adata, key_added=multiVI_CLUSTERS_KEY, resolution=res)
sc.pl.umap(adata, color = multiVI_CLUSTERS_KEY, legend_loc='on data', save="_leiden_fetalNephron_only.pdf", show=True)


# In[10]:


adata_UB = adata[adata.obs[multiVI_CLUSTERS_KEY].isin(["37","39"])]
adata_UB.obs['annoCellType'] = adata_UB.obs['type'].astype(str) + '_' + adata_UB.obs['cellType'].astype(str)
print(adata_UB.obs['annoCellType'].value_counts())
cell_types = adata_UB.obs['annoCellType'].to_list()


# In[11]:


to_remove = ["0","8","13","20","23","37","39","45","47"]

adata = adata[~adata.obs[multiVI_CLUSTERS_KEY].isin(to_remove)]
sc.pl.umap(adata, color = multiVI_CLUSTERS_KEY, legend_loc='on data', save="_leiden_fetalNephron_only_cleaned.pdf", show=True)

sc.pp.neighbors(adata, use_rep=multiVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=min_distUMAP)
sc.pl.umap(adata, color = "cellType", legend_loc='on data', save="_cellType_fetalNephron_only_cleaned.pdf", show=True)

adata.write_h5ad(filename = fetal_only_multiOme_nephron)

