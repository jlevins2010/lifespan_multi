#!/usr/bin/env python
# coding: utf-8

# ### Import libraries

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata
import scanpy as sc
import scvi
import matplotlib as mpl
import torch
import os.path
from IPython.display import display


# ### SCVI Settings

# In[2]:


scvi.__version__

scvi.settings.seed = 0
batch_size = 64
torch.set_float32_matmul_precision('high')


# ### Scanpy settings

# In[3]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450


# ### Sample Demographics

# In[4]:


input_file = "/vast/projects/chenyuli/nephrobase/levinsj/adata/HK3547_initialQC.h5ad"
output_file = "/vast/projects/chenyuli/nephrobase/levinsj/adataHK3547_postSOLO.h5ad"
if os.path.isfile(input_file):
    if not os.path.isfile(output_file):
        try:
                adata = sc.read_h5ad(input_file)
                print(adata)
                adata.var_names_make_unique()
                adata.X = adata.layers["counts"]
                scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key=None, labels_key=None) 
                vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb", dropout_rate=0.1)
                vae.view_anndata_setup(adata)
                vae.train(max_epochs = 1000, plan_kwargs={"lr":0.0005}, early_stopping = True, batch_size=batch_size) 
                solo = scvi.external.SOLO.from_scvi_model(vae)
                solo.train(batch_size=batch_size)
                df = solo.predict()
                df["prediction"] = solo.predict(soft= False)
                adata.obs["SOLO"] = df.prediction.astype(str)
                adata.obs = adata.obs.drop(columns=['_scvi_batch', '_scvi_labels'], errors='ignore') 
                adata = adata.copy() 
                adata.write(output_file) 
        except Exception as e:
                print("An unexpected error occurred:", e)


# In[5]:


input_file = "/vast/projects/chenyuli/nephrobase/levinsj/adata/HK3548_initialQC.h5ad"
output_file = "/vast/projects/chenyuli/nephrobase/levinsj/adataHK3548_postSOLO.h5ad"
if os.path.isfile(input_file):
    if not os.path.isfile(output_file):
        try:
                adata = sc.read_h5ad(input_file)
                print(adata)
                adata.var_names_make_unique()
                adata.X = adata.layers["counts"]
                scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key=None, labels_key=None) 
                vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb", dropout_rate=0.1)
                vae.view_anndata_setup(adata)
                vae.train(max_epochs = 1000, plan_kwargs={"lr":0.0005}, early_stopping = True, batch_size=batch_size) 
                solo = scvi.external.SOLO.from_scvi_model(vae)
                solo.train(batch_size=batch_size)
                df = solo.predict()
                df["prediction"] = solo.predict(soft= False)
                adata.obs["SOLO"] = df.prediction.astype(str)
                adata.write(output_file) 
        except Exception as e:
                print("An unexpected error occurred:", e)

