#!/usr/bin/env python
# coding: utf-8

# ### Import Libraries

# In[1]:


import numpy as np
import pandas as pd
import os
import scanpy as sc
from IPython.display import display
import matplotlib as mpl
import matplotlib.pyplot as plt


# ### Scanpy settings

# In[2]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2
res = 1.0
n_features = 250000


# ### SCVI settings

# In[3]:


SCVI_LATENT_KEY = "X_scvi"
SCVI_CLUSTERS_KEY = "clusters_scvi"


# In[4]:


input_file = "/vast/projects/chenyuli/nephrobase/levinsj/adata/postSCVI_snRNA_preclean1_human.h5ad"
outFile = "/vast/projects/chenyuli/nephrobase/levinsj/adata/postSCVI_snRNA_postclean1_human.h5ad"


# ### load input file

# In[5]:


adata = sc.read_h5ad(input_file)
print(adata.obs_names.duplicated().any())

sc.pl.umap(adata, color = "PrepType")
sc.pl.umap(adata, color = "study")
sc.pl.umap(adata, color = "SOLO")
sc.pl.umap(adata, color = "Amulet_doublet")

sc.tl.leiden(adata, key_added=SCVI_CLUSTERS_KEY, resolution=res)
sc.pl.umap(adata, color = SCVI_CLUSTERS_KEY, legend_loc='on data', frameon = False, legend_fontsize=10, legend_fontoutline=2)


# In[6]:


bad_clusters = ['41']


# ### annotate doublets

# In[7]:


df = adata.obs[SCVI_CLUSTERS_KEY].value_counts().to_frame()
df.index = adata.obs[SCVI_CLUSTERS_KEY].value_counts().index

df['sampleBreakdown'] = object
df['soloBreakdown'] = object
df['amuletBreakdown'] = object

samples = list(set(adata.obs["sample"].values.ravel().tolist()))

for i in df.index:
    counts_sample = []
    SOLO_sample = []
    amulet_sample = []
    for j in samples:
        counts_sample.append(adata.obs[(adata.obs['sample'] == j) & (adata.obs[SCVI_CLUSTERS_KEY] == i)].shape[0]/adata.obs[(adata.obs[SCVI_CLUSTERS_KEY] == i)].shape[0])
    SOLO_sample.append(adata.obs[(adata.obs['SOLO'] == "doublet") & (adata.obs[SCVI_CLUSTERS_KEY] == i)].shape[0]/adata.obs[(adata.obs[SCVI_CLUSTERS_KEY] == i)].shape[0])
    amulet_sample.append(adata.obs[(adata.obs['Amulet_doublet'] == "doublet") & (adata.obs[SCVI_CLUSTERS_KEY] == i)].shape[0]/adata.obs[(adata.obs[SCVI_CLUSTERS_KEY] == i)].shape[0])
    df.at[i,'sampleBreakdown'] = counts_sample
    df.at[i,'soloBreakdown'] = SOLO_sample
    df.at[i,'amuletBreakdown'] = amulet_sample

print(df)


# ### Make a bar graph showing samples per cluster

# In[8]:


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


# In[9]:


bad_clusters.append(['34','36','38','41'])


# ### Plot solo doublets by cluster

# In[10]:


width = 0.9

# Loop over each cell type and create a stacked bar plot of the "sample" category
for i, cell_type in enumerate(df.index):
    ind = np.arange(1)
    S1 = df.soloBreakdown[i]

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


# In[11]:


bad_clusters.append(['16','17','20','26','30','40','41'])


# ### Plot Amulet doublets by cluster

# In[12]:


width = 0.9

# Loop over each cell type and create a stacked bar plot of the "sample" category
for i, cell_type in enumerate(df.index):
    ind = np.arange(1)
    S1 = df.amuletBreakdown[i]

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


bad_clusters.append('26')


# In[14]:


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


# ### Remove bad clusters

# In[15]:


leiden = adata.obs[SCVI_CLUSTERS_KEY].unique()
leiden_mapping = {cluster: cluster not in bad_clusters for cluster in leiden}
print(leiden_mapping)


# In[16]:


adata.obs["clean"] = adata.obs[SCVI_CLUSTERS_KEY]
adata.obs["clean"] = adata.obs[SCVI_CLUSTERS_KEY].map(leiden_mapping).astype('category')

sc.pl.umap(adata, color = "clean")

print(adata.obs["clean"].value_counts())

adata = adata[adata.obs["clean"] == True]

print(adata)


# In[17]:


keep_str = ["sample",'Amulet_doublet','SOLO','PrepType','study','type',SCVI_CLUSTERS_KEY,'clean']
keep_numer = ['n_genes_by_counts', 'total_counts','pct_counts_mt']

keep = keep_str + keep_numer
adata.obs= adata.obs[keep]

# ensure proper data types
for i in keep_str:
    adata.obs[i]= adata.obs[i].astype('str')

for i in keep_numer:
    adata.obs[i]= adata.obs[i].astype('float')

adata.write_h5ad(outFile)
print(adata.obs["sample"].value_counts())

