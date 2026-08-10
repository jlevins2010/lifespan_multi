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
n_features = 250000

black_list_file = '/home/levinsj/Applications/hg38-blacklist.bed'


# ### SCVI settings

# In[3]:


print(scvi.__version__)

os.environ["CUDA_VISIBLE_DEVICES"] = "0"
scvi.settings.seed = 0
max_epochs = 40
batch_size = 64
max_epochs = 60
earlyStopping = True
learning_rate = 0.00001
scvi_layers = 2
scvi_latent = 30
fraction_expressed = 0.01

PEAKVI_LATENT_KEY = "X_peakvi"
PEAKVI_CLUSTERS_KEY = "clusters_peakvi"


# In[4]:


outputdir = "/home/levinsj/MultiOme/ATAC_data/Human/03_integrated_w_doublets/"
h5ad_inputdir = "/home/levinsj/MultiOme/ATAC_data/Human/02_doublets_called/"
adataSet = "/home/levinsj/MultiOme/ATAC_data/Human/03_integrated_w_doublets/preSCVI_w_doublets_data.h5ads"
outFile = "/home/levinsj/MultiOme/ATAC_data/Human/04_SCVI/postSCVI_ATAC_preclean1.h5ad"
peakVI_model = "/home/levinsj/MultiOme/scvi_Models/all_snATAC_with_doublets_human_preclean1/"


# ### select features and plot pre-integration with called doublets

# ### train model

# In[5]:


adata = sc.read_h5ad(outFile)


# ### will need to retrain when get more samples...
# model = scvi.model.PEAKVI.load("/home/levinsj/MultiOme/scvi_Models/peakVI_temp/", adata) 
# model.train(max_epochs = 40, plan_kwargs={"lr":learning_rate}, early_stopping = True, batch_size=batch_size) 
# 
# model.save(peakVI_model, overwrite=True)
# 

# In[6]:


model = scvi.model.PEAKVI.load(peakVI_model, adata)  # Load, providing the same AnnData

latent = model.get_latent_representation()
adata.obsm[PEAKVI_LATENT_KEY] = latent
latent.shape

sc.pp.neighbors(adata, use_rep=PEAKVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=min_distUMAP)
sc.pl.umap(adata, color = "sample")
sc.pl.umap(adata, color = "type")


# model = scvi.model.PEAKVI.load(peakVI_model, adata=adata)
# 

# ### plot ELBO curve

# In[7]:


ax=model.history['elbo_train'].plot()


# ### Get latent data

# In[8]:


sc.pp.neighbors(adata, use_rep=PEAKVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=min_distUMAP)
sc.pl.umap(adata, color = "PrepType")
sc.pl.umap(adata, color = "study")
sc.pl.umap(adata, color = "type")


# In[9]:


adata.write(outFile)

