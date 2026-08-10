#!/usr/bin/env python
# coding: utf-8

# ### Import Libraries

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata as ad
import scanpy as sc
import scvi
import matplotlib as mpl
import torch
import os
import os.path
from IPython.display import display


# ### SCVI import

# In[2]:


scvi.__version__

scvi.settings.seed = 0
batch_size = 64
max_epochs = 2000
earlyStopping = True
learning_rate = 0.0005
scvi_layers = 3
scvi_latent = 30

torch.set_float32_matmul_precision('high')

SCVI_LATENT_KEY = "X_scvi"
SCVI_CLUSTERS_KEY = "clusters_scvi"


# ### SCANPY Settings

# In[3]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2


# In[4]:


scVI_model_update = "/vast/projects/chenyuli/nephrobase/levinsj/scvi_models/all_snRNA_with_doublets2_human/"
input_file = "/vast/projects/chenyuli/nephrobase/levinsj/adata/postSCVI_snRNA_postclean1_human.h5ad"
output_file = "/vast/projects/chenyuli/nephrobase/levinsj/adata/postSCVI_snRNA_preclean2_human.h5ad"


# ### demographics

# In[5]:


adata_merge = sc.read_h5ad(input_file)


# ### SCVI integrate

# In[6]:


adata_merge.var_names_make_unique()

adata_merge.X = adata_merge.layers["counts"]
scvi.model.SCVI.setup_anndata(
    adata_merge,
    batch_key="sample",
    layer="counts", categorical_covariate_keys=["study", "PrepType"],
    continuous_covariate_keys=["total_counts","pct_counts_mt"])

model = scvi.model.SCVI(adata_merge, n_layers=scvi_layers, n_latent=scvi_latent, gene_likelihood="nb", dropout_rate=0.1)
model.view_anndata_setup(adata_merge)

model.train(max_epochs = max_epochs, plan_kwargs={"lr":learning_rate}, early_stopping = earlyStopping, batch_size=batch_size) 

model.save(scVI_model_update, overwrite=True)


# ### Plot ELBO

# In[7]:


ax=model.history['elbo_train'].plot()
pl=model.history['elbo_validation'].plot(ax=ax)


# ### plot integration

# In[8]:


latent = model.get_latent_representation()
adata_merge.obsm[SCVI_LATENT_KEY] = latent
latent.shape


# In[9]:


sc.pp.neighbors(adata_merge, use_rep=SCVI_LATENT_KEY)
sc.tl.umap(adata_merge, min_dist=min_distUMAP)
sc.pl.umap(adata_merge, color = "SOLO")
sc.pl.umap(adata_merge, color = "Amulet_doublet")
sc.pl.umap(adata_merge, color = "sample")
sc.pl.umap(adata_merge, color = "type")

adata_merge.write(output_file)

