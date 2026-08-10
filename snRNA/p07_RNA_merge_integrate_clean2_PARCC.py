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
res = 2.0
n_features = 250000


# ### SCVI settings

# In[3]:


SCVI_LATENT_KEY = "X_scvi"
SCVI_CLUSTERS_KEY = "clusters_scvi"


# ### Demographics

# In[4]:


input_file = "/vast/projects/chenyuli/nephrobase/levinsj/adata/postSCVI_snRNA_preclean2_human.h5ad"
output_file = "/vast/projects/chenyuli/nephrobase/levinsj/adata/postSCVI_snRNA_postclean2_human.h5ad"
barcodes_output_file = "/vast/projects/chenyuli/nephrobase/levinsj/adata/filtered_snRNA_Barcodes/postscVI_barcodes_human.csv"


# ### load input file

# In[5]:


adata = sc.read_h5ad(input_file)
sc.pl.umap(adata, color = "PrepType")
sc.pl.umap(adata, color = "study")
sc.pl.umap(adata, color = "SOLO")
sc.pl.umap(adata, color = "Amulet_doublet")

sc.tl.leiden(adata, key_added=SCVI_CLUSTERS_KEY, resolution=res)
sc.pl.umap(adata, color = SCVI_CLUSTERS_KEY, legend_loc='on data', frameon = False, legend_fontsize=10, legend_fontoutline=2)


# In[6]:


bad_clusters = []


# ### annotate doublets

# In[7]:


print(adata.obs_names[0:10])
print(adata.obs_names[-10:])


# In[8]:


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


# ### Make a bar graph showing samples per cluster

# In[9]:


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


# In[10]:


bad_clusters.append(['56'])


# ### Plot solo doublets by cluster

# In[11]:


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


# In[12]:


bad_clusters.append(["30","40","52"])


# ### Plot Amulet doublets by cluster

# In[13]:


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


# In[14]:


bad_clusters.append(['52'])


# In[15]:


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

# In[16]:


leiden = adata.obs[SCVI_CLUSTERS_KEY].unique()
leiden_mapping = {cluster: cluster not in bad_clusters for cluster in leiden}
print(leiden_mapping)


# In[17]:


adata.obs["clean"] = adata.obs[SCVI_CLUSTERS_KEY]
adata.obs["clean"] = adata.obs[SCVI_CLUSTERS_KEY].map(leiden_mapping).astype('category')

sc.pl.umap(adata, color = "clean")

print(adata.obs["clean"].value_counts())

adata = adata[adata.obs["clean"] == True]
adata = adata[adata.obs["SOLO"] != "doublet"]
adata = adata[adata.obs["Amulet_doublet"] != "doublet"]

print(adata)


# In[18]:


sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=min_distUMAP)
sc.tl.leiden(adata, key_added=SCVI_CLUSTERS_KEY, resolution=res)
sc.pl.umap(adata, color = SCVI_CLUSTERS_KEY, legend_loc='on data', frameon = False, legend_fontsize=10, legend_fontoutline=2)

sc.pl.umap(adata, color = "PrepType")
sc.pl.umap(adata, color = "study")
sc.pl.umap(adata, color = "SOLO")
sc.pl.umap(adata, color = "type")
sc.pl.umap(adata, color = "Amulet_doublet")


# ### Plot markers

# In[19]:


NPC_markers =  ["UNCX","SIX2","SIX1","CITED1"]
PT_markers  =  ["CUBN","APOE"]
iPT_markers  = ["HAVCR1","VCAM1"]
Int_markers =  ["JAG1","LHX1"]
PEC_markers =  ["CFH","WT1","PTPRO"]
DCT_markers =  ["WNK1","AQP3","SLC12A3","TFAP2A"]
LOH_markers =  ["UMOD"]
Podo_markers = ["MAFB","NPHS2","NPHS1"]
Stro_markers = ["COL3A1","COL1A1",'PDGFRA']
Endo_markers = ["EFGL7", "PLVAP","EMCN"]


def keep_matching_elements(my_list, adata_var_names):
  return [element for element in my_list if element in adata_var_names]

sc.pl.umap(adata, color = keep_matching_elements(NPC_markers, adata.var_names), color_map = "viridis_r")
sc.pl.umap(adata, color = keep_matching_elements(iPT_markers, adata.var_names), color_map = "viridis_r")
sc.pl.umap(adata, color = keep_matching_elements(PT_markers, adata.var_names), color_map = "viridis_r")
sc.pl.umap(adata, color = keep_matching_elements(Int_markers, adata.var_names), color_map = "viridis_r")
sc.pl.umap(adata, color = keep_matching_elements(PEC_markers, adata.var_names), color_map = "viridis_r")
sc.pl.umap(adata, color = keep_matching_elements(PEC_markers, adata.var_names), color_map = "viridis_r")
sc.pl.umap(adata, color = keep_matching_elements(DCT_markers, adata.var_names), color_map = "viridis_r")
sc.pl.umap(adata, color = keep_matching_elements(LOH_markers, adata.var_names), color_map = "viridis_r")
sc.pl.umap(adata, color = keep_matching_elements(Podo_markers, adata.var_names), color_map = "viridis_r")
sc.pl.umap(adata, color = keep_matching_elements(Stro_markers, adata.var_names), color_map = "viridis_r")
sc.pl.umap(adata, color = keep_matching_elements(Endo_markers, adata.var_names), color_map = "viridis_r")


# ### Clean dataType

# In[20]:


keep_str = ["sample",'Amulet_doublet','SOLO','PrepType','study','type',SCVI_CLUSTERS_KEY,'clean']
keep_numer = ['n_genes_by_counts', 'total_counts','pct_counts_mt']

keep = keep_str + keep_numer
adata.obs= adata.obs[keep]

# ensure proper data types
for i in keep_str:
    adata.obs[i]= adata.obs[i].astype('str')

for i in keep_numer:
    adata.obs[i]= adata.obs[i].astype('float')

adata.write_h5ad(output_file)
print(adata.obs["sample"].value_counts())

