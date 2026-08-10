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


demographics = pd.read_csv('/home/levinsj/MultiOme/Sample_demographics.csv') 
demographics.index = demographics["SampleName"]
display(demographics)


# In[5]:


input_file =  "/home/levinsj/MultiOme/CellBenderCorrected_snRNA/02_SCVI/postSCVI_snRNA_postclean2_human.h5ad"
output_file = "/home/levinsj/MultiOme/CellBenderCorrected_snRNA/02_SCVI/postSCVI_snRNA_annotated_human.h5ad"


# ### load input file

# In[6]:


adata = sc.read_h5ad(input_file)
sc.pl.umap(adata, color = "PrepType")
sc.pl.umap(adata, color = "study")
sc.pl.umap(adata, color = "SOLO")
sc.pl.umap(adata, color = "Amulet_doublet")
sc.pl.umap(adata, color = "type")


sc.tl.leiden(adata, key_added=SCVI_CLUSTERS_KEY, resolution=res)
sc.pl.umap(adata, color = SCVI_CLUSTERS_KEY, legend_loc='on data', frameon = False, legend_fontsize=10, legend_fontoutline=2)


# ### Plot markers

# In[7]:


NPC_markers =  ["UNCX","SIX2","SIX1","CITED1"]
PT_markers  =  ["CUBN","APOE"]
Int_markers =  ["JAG1","LHX1"]
PEC_markers =  ["CFH","WT1","PTPRO"]
DCT_markers =  ["WNK1","AQP3","SLC12A3","TFAP2A"]
LOH_markers =  ["UMOD","AQP1","SLC12A1"]
Podo_markers = ["MAFB","NPHS2","NPHS1"]
Stro_markers = ["COL3A1","COL1A1",'PDGFRA']
Endo_markers = ["EFGL7", "PLVAP","EMCN"]
CNT_markers = ["TRPV5","CALB1","VDR","KL"]
IC_markers = ["SLC4A1",'Adgrf5','SLC26A7','SLC26A4','ATP6V1B1']
UB_markers = ["RET",'GATA3']


def keep_matching_elements(my_list, adata_var_names):
  return [element for element in my_list if element in adata_var_names]


sc.pl.umap(adata, color = keep_matching_elements(Stro_markers, adata.var_names), color_map = "viridis_r")
sc.pl.umap(adata, color = keep_matching_elements(Endo_markers, adata.var_names), color_map = "viridis_r")


# In[8]:


sc.pl.umap(adata, color = keep_matching_elements(NPC_markers, adata.var_names), color_map = "viridis_r")
sc.pl.umap(adata, color = keep_matching_elements(Int_markers, adata.var_names), color_map = "viridis_r")


# In[9]:


sc.pl.umap(adata, color = keep_matching_elements(Podo_markers, adata.var_names), color_map = "viridis_r")
sc.pl.umap(adata, color = keep_matching_elements(PT_markers, adata.var_names), color_map = "viridis_r")
sc.pl.umap(adata, color = keep_matching_elements(PEC_markers, adata.var_names), color_map = "viridis_r")


# In[10]:


sc.pl.umap(adata, color = keep_matching_elements(DCT_markers, adata.var_names), color_map = "viridis_r")
sc.pl.umap(adata, color = keep_matching_elements(CNT_markers, adata.var_names), color_map = "viridis_r")
sc.pl.umap(adata, color = keep_matching_elements(LOH_markers, adata.var_names), color_map = "viridis_r")


# In[11]:


sc.pl.umap(adata, color = keep_matching_elements(IC_markers, adata.var_names), color_map = "viridis_r")


# In[12]:


adata2 = adata.copy()
adata2.X = adata2.layers["counts"] # set X-layer prior to subset
sc.pp.normalize_total(adata2, target_sum=1e4)
sc.pp.log1p(adata2)
sc.tl.rank_genes_groups(adata2, SCVI_CLUSTERS_KEY, method='wilcoxon')
sc.pl.rank_genes_groups(adata2, n_genes=25, sharey=False)


# In[13]:


adata.obs["cellType"] = adata.obs[SCVI_CLUSTERS_KEY]
adata = adata[~adata.obs[SCVI_CLUSTERS_KEY].isin(['10','26','56'])] ### remove neurons, neurons and *** respectively

cell_identities = {'0':'LOH','1':'DCT_CNT','2':'PT','3':'LOH','4':'DCT_CNT','5':'DCT_CNT','6':'LOH',
                   '7':'Stroma','8':'Stroma','9':'PT','11':'Endothelium','12':'PT','13':'Podocyte','14':'Endothelium',
                   '15':'PT', '16':'LOH','17':'NPC','18':'NPC', '19':'OT','20':'PT','21':'IC','22':'PT',
                   '23':'LOH','24':'PT','25':'PT','27':'PEC', '28':'Stroma','29':'Stroma','30':'LOH',
                   '31':'Int','32':'IC', '33': 'PT', '34':'Int','35':'Immune Cells', '36':'Int','37':'Immune Cells',
                   '38':'DCT_CNT','39':'LOH', '40': 'Endothelium', '41':'Stroma','42':'Stroma','43':'Int','44':'DCT_CNT',
                   '45':'DCT_CNT','46':'DCT_CNT','47':'LOH','48':'Immune Cells', '49':'LOH','50':'DCT_CNT','51':'Stroma','52':'Endothelium',
                   '53':'Stroma', '54':'LOH', '55':'Endothelium', '57': 'NPC'}
adata.obs["cellType_RNA"] = adata.obs[SCVI_CLUSTERS_KEY].map(cell_identities).astype('category')
sc.pl.umap(adata, color = "cellType_RNA")


# In[14]:


shared_col = demographics.columns.get_loc('snRNA_finalCells')

for index, item in enumerate(demographics['sampleRawH5Matricies_snRNA']):
    sampleName = demographics['SampleName'][index]
    if demographics['SampleName'][index] in adata.obs["sample"].values:
        adata_subset = adata[adata.obs["sample"] == sampleName]
        demographics.iloc[index, shared_col] = adata_subset.n_obs
    else:
        demographics.iloc[index, shared_col] = 0

demographics.to_csv('/home/levinsj/MultiOme/Sample_demographics.csv', index=False)


# In[15]:


sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=min_distUMAP)
sc.tl.leiden(adata, key_added=SCVI_CLUSTERS_KEY, resolution=res)
sc.pl.umap(adata, color = SCVI_CLUSTERS_KEY, legend_loc='on data', frameon = False, legend_fontsize=10, legend_fontoutline=2)


# In[16]:


adata.write(output_file)

