#!/usr/bin/env python
# coding: utf-8

# ### Import libraries

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata
import scanpy as sc
import matplotlib as mpl
import os.path
from IPython.display import display


# ### Scanpy Settings

# In[2]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450


# ### demographics

# In[3]:


demographics = pd.read_csv('/home/levinsj/MultiOme/Sample_demographics.csv') 
demographics.index = demographics["SampleName"]
display(demographics)


# ### mark snRNAseq Doublets in snATAC data

# In[4]:


for index, item in enumerate(demographics['sampleRawH5Matricies_snRNA']):
    input_file = "/home/levinsj/MultiOme/CellBenderCorrected_snRNA/03_post_SOLO/"+demographics['SampleName'][index]+"_postSOLO.h5ad"
    output_file = "/home/levinsj/MultiOme/CellBenderCorrected_snRNA/04_postAmulet/"+demographics['SampleName'][index]+"_postSOLO.h5ad"
    amuletDoublet = "/home/levinsj/MultiOme/Barcodes/amulet_doublets/" +demographics['SampleName'][index]+ "_MultipletCellIds_01.txt"
    if os.path.isfile(input_file):
        if demographics['PrepType'][index] == "MultiOme":
            if os.path.isfile(amuletDoublet):
                # open the amulet doublet file and return barcodes that are doublets
                print("Checking amulet annotations for " + demographics['SampleName'][index])
                adata = sc.read_h5ad(input_file)
                adata.obs["PrepType"] = "multiOme"
                csv_df = pd.read_csv(amuletDoublet)
                csv_index = set(csv_df.iloc[:, 0])
                adata.obs['Amulet_doublet'] = adata.obs.index.isin(csv_index).astype(str)
                adata.obs['Amulet_doublet'] = adata.obs['Amulet_doublet'].replace({"True": 'doublet', 'False': 'singlet'})
                adata.write(output_file)
                print("Perportion singlets: " + str(adata.obs[(adata.obs['Amulet_doublet'] == "singlet") & (adata.obs['SOLO'] == 'singlet')].shape[0]/adata.n_obs))
                print("Amulet doublet, SOLO singlet: " + str(adata.obs[(adata.obs['Amulet_doublet'] == "doublet") & (adata.obs['SOLO'] == 'singlet')].shape[0]/adata.n_obs))
                print("Amulet singlet, SOLO doublet: " + str(adata.obs[(adata.obs['Amulet_doublet'] == "singlet") & (adata.obs['SOLO'] == 'doublet')].shape[0]/adata.n_obs))
                print("Amulet doublet, SOLO doublet: " + str(adata.obs[(adata.obs['Amulet_doublet'] == "doublet") & (adata.obs['SOLO'] == 'doublet')].shape[0]/adata.n_obs))
            else:
                print("No amulet file for " + demographics['SampleName'][index])
        else:
            print(demographics['SampleName'][index] + " is not a multiOme! No amulet doublets to add.")
            adata = sc.read_h5ad(input_file)
            adata.obs["Amulet_doublet"] = np.nan
            adata.obs["PrepType"] = "snRNA"
            adata.write(output_file)
    else:
        print(demographics['SampleName'][index] + " does not have a post SOLO file!")

