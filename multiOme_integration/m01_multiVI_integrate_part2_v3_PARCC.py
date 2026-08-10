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
import torch
import matplotlib as mpl
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
learning_rate = 0.0001
scvi_layers = 3
scvi_latent = 30
fraction_expressed = 0.01

torch.set_float32_matmul_precision('high')

multiVI_LATENT_KEY = "X_multivi"
multiVI_CLUSTERS_KEY = "clusters_multivi"


# ### SCANPY Settings

# In[3]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2


# ### demographics

# In[4]:


demographics = pd.read_csv('/vast/projects/chenyuli/nephrobase/levinsj/input_csv/Sample_demographics.csv') 
demographics.index = demographics["SampleName"]
display(demographics)


# ### Get shared barcodes

# In[5]:


snRNAinput = "/vast/projects/chenyuli/nephrobase/levinsj/adata/postSCVI_snRNA_annotated_human.h5ad"
snATACinput = "/vast/projects/chenyuli/nephrobase/levinsj/adata/postSCVI_ATAC_peaks_called.h5ad"

preSCVIfile = "/vast/projects/chenyuli/nephrobase/levinsj/adata/merged_preMultiVI_human.h5ad"
multiVImodel = "/vast/projects/chenyuli/nephrobase/levinsj/scvi_models/multiVI_pretrained_integration1/"
multiVImodel_final = "/vast/projects/chenyuli/nephrobase/levinsj/scvi_models/multiVI_posttrained_integration1/"
preSCVIfile_filtered = "/vast/projects/chenyuli/nephrobase/levinsj/adata/merged_preMultiVI_human_filtered.h5ad"

postSCVIfile = "/vast/projects/chenyuli/nephrobase/levinsj/adata/merged_postMultiVI_human.h5ad"


# In[6]:


torch.serialization.add_safe_globals([np.core.multiarray._reconstruct])


# Combine the multiOme into one object

# In[7]:


adata_mvi = sc.read_h5ad(preSCVIfile)
sc.pp.filter_genes(adata_mvi, min_cells=int(adata_mvi.shape[0] * fraction_expressed))
adata_mvi = adata_mvi.copy() 

adata_mvi.write_h5ad(filename = preSCVIfile_filtered)
print(adata_mvi.shape)


# ### SCVI integrate

# In[8]:


scvi.model.MULTIVI.setup_anndata(
    adata_mvi, 
    batch_key="modality",
    categorical_covariate_keys=["sample", "PrepType", "type"],
    continuous_covariate_keys=["total_counts", "pct_counts_mt"]
)

print(adata_mvi.obs["modality"].value_counts())

# 2. Immediately build the model using the exact same object
model = scvi.model.MULTIVI(
    adata_mvi,
    n_genes=(adata_mvi.var["modality"] == "1").sum(),
    n_regions=(adata_mvi.var["modality"] == "0").sum(),
)

model.view_anndata_setup()
model.save(multiVImodel, overwrite=True)


# In[9]:


adata_mvi = sc.read_h5ad(preSCVIfile_filtered)
print(adata_mvi)


# In[10]:


model.train(max_epochs=500, plan_kwargs={"lr":learning_rate}, early_stopping = earlyStopping, batch_size=batch_size)

model.save(multiVImodel_final, overwrite=True)

latent = model.get_latent_representation()
adata_mvi.obsm[multiVI_LATENT_KEY] = latent
latent.shape

sc.pp.neighbors(adata_mvi, use_rep=multiVI_LATENT_KEY)
sc.tl.umap(adata_mvi, min_dist=min_distUMAP)
sc.pl.umap(adata_mvi, color = "sample")
sc.pl.umap(adata_mvi, color = "modality")
sc.pl.umap(adata_mvi, color = "type")


# ### Plot ELBO

# In[11]:


ax=model.history['elbo_train'].plot()
pl=model.history['elbo_validation'].plot(ax=ax)


# ### plot integration

# In[12]:


adata_mvi.write(postSCVIfile)

