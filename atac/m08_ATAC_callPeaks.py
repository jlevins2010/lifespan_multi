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


# ### SNAPATAC settings

# In[2]:


snap.__version__

sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2
n_features = 250000

black_list_file = '/home/levinsj/Applications/hg38-blacklist.bed'


# ### Scanpy settings

# In[3]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2
res = 1.0
n_features = 250000


# ### SCVI settings

# In[4]:


PEAKVI_LATENT_KEY = "X_peakvi"
PEAKVI_CLUSTERS_KEY = "clusters_peakvi"


# ### Demographics

# In[5]:


demographics = pd.read_csv('/home/levinsj/MultiOme/Sample_demographics.csv') 
demographics.index = demographics["SampleName"]
display(demographics)


# In[6]:


input_anndataset = "/home/levinsj/MultiOme/ATAC_data/Human/03_integrated_w_doublets/preSCVI_w_doublets_data.h5ads"
output_anndataset = "/home/levinsj/MultiOme/ATAC_data/Human/03_integrated_w_doublets/postSCVI_final_data.h5ads"
input_filtered_adata = "/home/levinsj/MultiOme/ATAC_data/Human/04_SCVI/postSCVI_ATAC_postclean2.h5ad"
output_file = "/home/levinsj/MultiOme/ATAC_data/Human/04_SCVI/postSCVI_ATAC_peaks_called.h5ad"
filtered_barcodes = "/home/levinsj/MultiOme/Barcodes/snapATAC_Barcodes/postPeakVI_barcodes_human.csv"
output_gene_accessibility_adata = "/home/levinsj/MultiOme/ATAC_data/Human/04_SCVI/postSCVI_geneAccessibility.h5ad"


# ### load input file

# In[7]:


# get barcodes from cells that passed QC and filtering
snapBarcodes = pd.read_csv(filtered_barcodes, header=0) 
filtered_indecies = snapBarcodes.iloc[:, 0].tolist() 
corresponding_clusters = snapBarcodes.iloc[:, 1].tolist() 

print(len(filtered_indecies))
print(len(corresponding_clusters))
print(filtered_indecies[0:10])


# In[8]:


adataset = snap.read_dataset(filename=input_anndataset)
print(adataset)


# In[9]:


### Subset to just filtred indecies
obs_names = pd.Series(adataset.obs_names)  
keepidx = obs_names.isin(filtered_indecies)
print(keepidx.value_counts())

peakVI_clusters = [str(x) for x in corresponding_clusters]


# In[10]:


filtered_anndata = adataset.subset(obs_indices = keepidx, out = output_anndataset)
filtered_anndata[0].obs[PEAKVI_CLUSTERS_KEY] = peakVI_clusters
print(filtered_anndata)
adataset.close()


# In[11]:


print(filtered_anndata[0].isbacked)
print(filtered_anndata[0])


# In[12]:


try:
    print(filtered_anndata[0].obs[PEAKVI_CLUSTERS_KEY].value_counts())

except Exception as e:
    print(f"An error occurred: {e}")


# In[13]:


try:
    gene_matrix = snap.pp.make_gene_matrix(filtered_anndata[0], gene_anno=snap.genome.hg38)
    gene_matrix.write_h5ad(output_gene_accessibility_adata)
    print(gene_matrix)
except Exception as e:
    print(f"An error occurred: {e}")


# In[14]:


try:
    gene_matrix = snap.pp.make_gene_matrix(filtered_anndata, gene_anno=snap.genome.hg38)
    gene_matrix.write_h5ad(output_gene_accessibility_adata)
    print(gene_matrix)
except Exception as e:
    print(f"An error occurred: {e}")


# try:
#     snap.tl.macs3(filtered_anndata[0], groupby= PEAKVI_CLUSTERS_KEY)
#     peaks = snap.tl.merge_peaks(filtered_anndata[0].uns['macs3'], snap.genome.hg38)
#     peak_mat = snap.pp.make_peak_matrix(filtered_anndata[0], use_rep=peaks['Peaks'])
#     print(peak_mat)
#     peak_mat.write_h5ad(output_file)
# except Exception as e:
#     print(f"An error occurred: {e}")
# 

# In[15]:


filtered_anndata[0].close()

