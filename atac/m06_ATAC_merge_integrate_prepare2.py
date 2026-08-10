#!/usr/bin/env python
# coding: utf-8

# ### Import Libraries

# In[1]:


import snapatac2 as snap
import numpy as np
import pandas as pd
import os
import magic
import scanpy as sc
import shutil
import scvi
import torch
from glob import glob
from tqdm import tqdm
from IPython.display import display
import matplotlib as mpl


# ### Scanpy and SNAPATAC2 settings

# In[2]:


snap.__version__

sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2
res = 0.2
n_features = 250000

black_list_file = '/home/levinsj/Applications/hg38-blacklist.bed'


# ### SCVI settings

# In[3]:


print(scvi.__version__)

os.environ["CUDA_VISIBLE_DEVICES"] = "0"
scvi.settings.seed = 0
max_epochs = 2000
batch_size = 64
max_epochs = 60
earlyStopping = True
learning_rate = 0.0005
scvi_layers = 2
scvi_latent = 30

PEAKVI_LATENT_KEY = "X_peakvi"
PEAKVI_CLUSTERS_KEY = "clusters_peakvi"


# ### Demographics

# In[4]:


demographics = pd.read_csv('/home/levinsj/MultiOme/Sample_demographics.csv') 
demographics.index = demographics["SampleName"]
display(demographics)


# In[5]:


input_file = "/home/levinsj/MultiOme/ATAC_data/Human/04_SCVI/postSCVI_ATAC_postclean1.h5ad"
peakVI_model = "/home/levinsj/MultiOme/scvi_Models/all_snATAC_with_doublets_postfirstclean_human/"
output_file = "/home/levinsj/MultiOme/ATAC_data/Human/04_SCVI/post_SCVI_clean2.h5ad"


# In[6]:


adata = sc.read_h5ad(input_file)
print(adata)


# ### train model

# In[7]:


print(adata)
print(adata.obs["PrepType"].value_counts())
print(adata.obs["study"].value_counts())
print(adata.obs["type"].value_counts())


# In[8]:


scvi.model.PEAKVI.setup_anndata(
    adata,
    batch_key="sample", categorical_covariate_keys=["study", "PrepType"])

model = scvi.model.PEAKVI(adata, n_latent = scvi_latent, n_layers_encoder = scvi_layers, n_layers_decoder= scvi_layers,  dropout_rate=0.1)
model.view_anndata_setup(adata)


# In[9]:


model.train(max_epochs = max_epochs, plan_kwargs={"lr":0.00001}, early_stopping = True, batch_size=batch_size) 

model.save(peakVI_model, overwrite=True)


# ### plot ELBO curve

# In[10]:


ax=model.history['elbo_train'].plot()


# ### Get latent data

# In[11]:


latent = model.get_latent_representation()
adata.obsm[PEAKVI_LATENT_KEY] = latent
latent.shape


# In[12]:


sc.pp.neighbors(adata, use_rep=PEAKVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=min_distUMAP)
sc.pl.umap(adata, color = "PrepType")
sc.pl.umap(adata, color = "study")
sc.pl.umap(adata, color = "type")


# In[13]:


adata.write(output_file)

