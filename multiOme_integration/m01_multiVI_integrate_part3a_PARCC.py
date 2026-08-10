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

# In[4]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2


# ### demographics

# In[5]:


demographics = pd.read_csv('/vast/projects/chenyuli/nephrobase/levinsj/input_csv/Sample_demographics.csv') 
demographics.index = demographics["SampleName"]
display(demographics)


# In[6]:


outputFile = "/vast/projects/chenyuli/nephrobase/levinsj/adata/merged_postMultiVI_V3_human_filtered.h5ad"
postSCVIfile = "/vast/projects/chenyuli/nephrobase/levinsj/adata/merged_postMultiVI_human.h5ad"


# In[7]:


adata = sc.read_h5ad(postSCVIfile)
bad_clusters = []


# In[8]:


sc.tl.leiden(adata, key_added=multiVI_CLUSTERS_KEY, resolution=res)
sc.pl.umap(adata, color = multiVI_CLUSTERS_KEY, legend_loc='on data', frameon = False, legend_fontsize=10, legend_fontoutline=2)


# In[9]:


df = adata.obs[multiVI_CLUSTERS_KEY].value_counts().to_frame()
df.index = adata.obs[multiVI_CLUSTERS_KEY].value_counts().index

df['sampleBreakdown'] = object
df['PrepType'] = object
df['modality'] = object

samples = list(set(adata.obs["sample"].values.ravel().tolist()))

for i in df.index:
    counts_sample = []
    prep_sample = []
    modality_sample = []
    for j in samples:
        counts_sample.append(adata.obs[(adata.obs['sample'] == j) & (adata.obs[multiVI_CLUSTERS_KEY] == i)].shape[0]/adata.obs[(adata.obs[multiVI_CLUSTERS_KEY] == i)].shape[0])
    prep_sample.append(adata.obs[(adata.obs['PrepType'] == "multiOme") & (adata.obs[multiVI_CLUSTERS_KEY] == i)].shape[0]/adata.obs[(adata.obs[multiVI_CLUSTERS_KEY] == i)].shape[0])
    modality_sample.append(adata.obs[(adata.obs['modality'] == "paired") & (adata.obs[multiVI_CLUSTERS_KEY] == i)].shape[0]/adata.obs[(adata.obs[multiVI_CLUSTERS_KEY] == i)].shape[0])
    df.at[i,'sampleBreakdown'] = counts_sample
    df.at[i,'PrepType'] = prep_sample
    df.at[i,'modality'] = modality_sample


# In[10]:


width = 0.9

# Loop over each cell type and create a stacked bar plot of the "sample" category
for i, cell_type in enumerate(df.index):
    ind = np.arange(1)
    S1 = df.sampleBreakdown[i]

    # generate a palette with all samples (scramble colors)
    num_colors = len(samples)
    cmap = plt.cm.get_cmap('hsv', num_colors)
    palette_indices = np.random.permutation(num_colors)
    palette = [cmap(i) for i in palette_indices]

    bottom = 0
    for j in range(len(S1)):
          plt.bar(i, S1[j],bottom = bottom, width = width, color = palette[j])
          bottom = bottom + S1[j]
    # Show the plot
plt.tight_layout()

plt.show()


# In[11]:


bad_clusters.append(['56'])


# In[12]:


width = 0.9

# Loop over each cell type and create a stacked bar plot of the "sample" category
for i, cell_type in enumerate(df.index):
    ind = np.arange(1)
    S1 = df.PrepType[i]

    # generate a palette with all samples (scramble colors)
    num_colors = 2
    cmap = plt.cm.get_cmap('hsv', num_colors)
    palette = [cmap(i) for i in palette_indices]

    bottom = 0
    for j in range(len(S1)):
          plt.bar(i, S1[j],bottom = bottom, width = width, color = palette[j])
          bottom = bottom + S1[j]
    # Show the plot
plt.tight_layout()

plt.show()


# In[13]:


bad_clusters.append([54])


# In[14]:


width = 0.9

# Loop over each cell type and create a stacked bar plot of the "sample" category
for i, cell_type in enumerate(df.index):
    ind = np.arange(1)
    S1 = df.modality[i]

    # generate a palette with all samples (scramble colors)
    num_colors = 2
    cmap = plt.cm.get_cmap('hsv', num_colors)
    palette = [cmap(i) for i in palette_indices]

    bottom = 0
    for j in range(len(S1)):
          plt.bar(i, S1[j],bottom = bottom, width = width, color = palette[j])
          bottom = bottom + S1[j]
    # Show the plot
plt.tight_layout()

plt.show()


# In[15]:


bad_clusters.append(['0','1','2','16','38','48','49','52','53','54'])


# In[16]:


def flatten(l):
  flat_list = []
  for sublist in l:
    if isinstance(sublist, list):
      for item in sublist:
        flat_list.append(item)
    else:
      flat_list.append(sublist)
  return flat_list

bad_clusters = (flatten(bad_clusters))
print(bad_clusters)


# In[17]:


leiden = adata.obs[multiVI_CLUSTERS_KEY].unique()
leiden_mapping = {cluster: cluster not in bad_clusters for cluster in leiden}
print(leiden_mapping)


# In[18]:


adata.obs["clean"] = adata.obs[multiVI_CLUSTERS_KEY]
adata.obs["clean"] = adata.obs[multiVI_CLUSTERS_KEY].map(leiden_mapping).astype('category')

sc.pl.umap(adata, color = "clean")
sc.pl.umap(adata, color = "type")
sc.pl.umap(adata, color = "modality")

print(adata.obs["clean"].value_counts())

adata = adata[adata.obs["clean"] == True]

sc.pl.umap(adata, color = "clean")
sc.pl.umap(adata, color = "type")
sc.pl.umap(adata, color = "modality")


# In[19]:


sc.pp.neighbors(adata, use_rep=multiVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=min_distUMAP)
sc.tl.leiden(adata, key_added=multiVI_CLUSTERS_KEY, resolution=res)
sc.pl.umap(adata, color = "sample")
sc.pl.umap(adata, color = "modality")
sc.pl.umap(adata, color = "type")
sc.pl.umap(adata, color = multiVI_CLUSTERS_KEY, legend_loc='on data', frameon = False, legend_fontsize=10, legend_fontoutline=2)


# In[20]:


adata = adata[~adata.obs[multiVI_CLUSTERS_KEY].isin(['9','17','37','41'])]
sc.pp.neighbors(adata, use_rep=multiVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=min_distUMAP)
sc.tl.leiden(adata, key_added=multiVI_CLUSTERS_KEY, resolution=res)
sc.pl.umap(adata, color = "sample")
sc.pl.umap(adata, color = "modality")
sc.pl.umap(adata, color = "type")
sc.pl.umap(adata, color = multiVI_CLUSTERS_KEY, legend_loc='on data', frameon = False, legend_fontsize=10, legend_fontoutline=2)


# In[21]:


keep_str = ["sample",'Amulet_doublet','SOLO','PrepType','study','type',multiVI_CLUSTERS_KEY,'clean']
keep_numer = ['n_genes_by_counts', 'total_counts','pct_counts_mt']

keep = keep_str + keep_numer
adata.obs= adata.obs[keep]

# ensure proper data types
for i in keep_str:
    adata.obs[i]= adata.obs[i].astype('str')

for i in keep_numer:
    adata.obs[i]= adata.obs[i].astype('float')

adata.write_h5ad(outputFile)
print(adata)

