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
batch_size = 32
max_epochs = 60
earlyStopping = True
learning_rate = 0.00001
scvi_layers = 2
scvi_latent = 30
fraction_expressed = 0.01

PEAKVI_LATENT_KEY = "X_peakvi"
PEAKVI_CLUSTERS_KEY = "clusters_peakvi"


# ### Demographics

# In[4]:


demographics = pd.read_csv('/home/levinsj/MultiOme/Sample_demographics.csv') 
demographics.index = demographics["SampleName"]
display(demographics)


# In[5]:


outputdir = "/home/levinsj/MultiOme/ATAC_data/Human/03_integrated_w_doublets/"
h5ad_inputdir = "/home/levinsj/MultiOme/ATAC_data/Human/02_doublets_called/"
adataSet = "/home/levinsj/MultiOme/ATAC_data/Human/03_integrated_w_doublets/preSCVI_w_doublets_data.h5ads"
outFile = "/home/levinsj/MultiOme/ATAC_data/Human/04_SCVI/postSCVI_ATAC_preclean1.h5ad"
peakVI_model = "/home/levinsj/MultiOme/scvi_Models/all_snATAC_with_doublets_human_10epochs/"


# ### make a directory with sample data for futher processing as an integrated object

# In[6]:


os.makedirs(outputdir, exist_ok=True)

adatas_to_transfer = glob(os.path.join(h5ad_inputdir, '*.h5ad'), recursive=False)
for fl in adatas_to_transfer:
    shutil.copy(fl, outputdir)


# In[7]:


adatas_to_read = glob(os.path.join(outputdir, '*.h5ad'), recursive=False)

print(adatas_to_read)


# In[8]:


adatas = []


for fl in adatas_to_read:
    adata = snap.read(filename=fl)
    adatas.append(adata)


# In[9]:


adataset = snap.AnnDataSet(
    adatas=[(f.filename.split('/')[-1].split('.h5ad')[0], f) for f in adatas],
    filename= adataSet
)

print(adataset)


# ### select features and plot pre-integration with called doublets

# In[10]:


snap.pp.select_features(adataset, n_features= n_features, blacklist=black_list_file)

unique_cell_ids = [sa + '_' + bc for sa, bc in zip(adataset.obs['sample'], adataset.obs_names)]
adataset.obs_names = unique_cell_ids
assert adataset.n_obs == np.unique(adataset.obs_names).size


# In[11]:


adataset.obs["Amulet_doublet"] = adataset.adatas.obs['Amulet_doublet']
adataset.obs["SOLO"] = adataset.adatas.obs['SOLO']


# ### train model

# In[12]:


adata = adataset.to_adata()
adata

adataset.close()


# In[13]:


demo_dict = dict(zip(demographics['SampleName'], demographics['PrepType']))
adata.obs["PrepType"] = adata.obs["sample"]
adata.obs["PrepType"] = adata.obs['sample'].map(demo_dict).astype('category')

demo_dict = dict(zip(demographics['SampleName'], demographics['Site']))
adata.obs["study"] = adata.obs["sample"]
adata.obs["study"] = adata.obs['sample'].map(demo_dict).astype('category')

demo_dict = dict(zip(demographics['SampleName'], demographics['Type']))
adata.obs["type"] = adata.obs["sample"]
adata.obs["type"] = adata.obs['sample'].map(demo_dict).astype('category')

print(adata)
print(adata.obs["PrepType"].value_counts())
print(adata.obs["study"].value_counts())
print(adata.obs["type"].value_counts())


# In[15]:


sc.pp.filter_genes(adata, min_cells=int(adata.shape[0] * fraction_expressed))

adata.write(outFile)

scvi.model.PEAKVI.setup_anndata(
    adata,
    batch_key="sample", categorical_covariate_keys=["study", "PrepType"])

model = scvi.model.PEAKVI(adata, n_latent = scvi_latent, n_layers_encoder = scvi_layers, n_layers_decoder= scvi_layers,  dropout_rate=0.1)
model.view_anndata_setup(adata)


# In[16]:


model.train(max_epochs = 10, plan_kwargs={"lr":learning_rate}, early_stopping = True, batch_size=batch_size) 
model.save("/home/levinsj/MultiOme/scvi_Models/peakVI_temp/", overwrite=True)

model = scvi.model.PEAKVI.load("/home/levinsj/MultiOme/scvi_Models/peakVI_temp/", adata) 
model.train(max_epochs = 10, plan_kwargs={"lr":learning_rate}, early_stopping = True, batch_size=batch_size) 
model.save("/home/levinsj/MultiOme/scvi_Models/peakVI_temp/", overwrite=True)

model = scvi.model.PEAKVI.load("/home/levinsj/MultiOme/scvi_Models/peakVI_temp/", adata) 
model.train(max_epochs = 10, plan_kwargs={"lr":learning_rate}, early_stopping = True, batch_size=batch_size) 
model.save("/home/levinsj/MultiOme/scvi_Models/peakVI_temp/", overwrite=True)

model = scvi.model.PEAKVI.load("/home/levinsj/MultiOme/scvi_Models/peakVI_temp/", adata) 
model.train(max_epochs = 10, plan_kwargs={"lr":learning_rate}, early_stopping = True, batch_size=batch_size) 
model.save("/home/levinsj/MultiOme/scvi_Models/peakVI_temp/", overwrite=True)

model = scvi.model.PEAKVI.load("/home/levinsj/MultiOme/scvi_Models/peakVI_temp/", adata) 
model.train(max_epochs = 10, plan_kwargs={"lr":learning_rate}, early_stopping = True, batch_size=batch_size) 
model.save("/home/levinsj/MultiOme/scvi_Models/peakVI_temp/", overwrite=True)

