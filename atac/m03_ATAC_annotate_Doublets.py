#!/usr/bin/env python
# coding: utf-8

# ### Import libraries

# In[1]:


import snapatac2 as snap
import numpy as np
import pandas as pd
import os
import magic
import scanpy as sc
from glob import glob
from tqdm import tqdm
from IPython.display import display
import matplotlib as mpl


# Need to run amulet first before using this script. To run amulet you need to
# 
# 
# (1) run run_amulet.sh python script

# ### Scanpy and SNAPATAC2 settings

# In[2]:


snap.__version__

sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450


# ### Demographics

# In[3]:


demographics = pd.read_csv('/home/levinsj/MultiOme/Sample_demographics.csv') 
demographics.index = demographics["SampleName"]
display(demographics)


# ### filter the adata object for QC and export barcodes of all cells called

# In[4]:


for index, item in enumerate(demographics['fragment_files_snATAC']):
    input_file = '/home/levinsj/MultiOme/ATAC_data/Human/01_filtered/' + demographics['SampleName'][index]+ ".h5ad"
    output_file = '/home/levinsj/MultiOme/ATAC_data/Human/02_doublets_called/' + demographics['SampleName'][index]+ ".h5ad"
    amuletDoublet = "/home/levinsj/MultiOme/Barcodes/amulet_doublets/" +demographics['SampleName'][index]+ "_MultipletCellIds_01.txt"
    soloDoublet = "/home/levinsj/MultiOme/Barcodes/SOLO_doublets/" +demographics['SampleName'][index]+ "_soloBarCodes.csv"
    if os.path.isfile(input_file):
        if os.path.isfile(amuletDoublet):
            print("Checking amulet annotations for " + demographics['SampleName'][index])
            adata_input = snap.read(filename=input_file)
            adata = adata_input.copy(output_file)
            csv_df = pd.read_csv(amuletDoublet)
            csv_index = list(csv_df.iloc[:, 0])
            in_csv_dict = {index: index in csv_index for index in adata.obs_names}
            adata.obs['Amulet_doublet'] = [str(in_csv_dict[name]) for name in adata.obs_names]
            adata.obs['Amulet_doublet'] = adata.obs['Amulet_doublet'].replace({"True": 'doublet', 'False': 'singlet'})
            adata.obs["study"] = [demographics['Site'][index]]  * len(adata.obs_names)
            adata.obs["sample"] = [demographics['SampleName'][index]]  * len(adata.obs_names)
            if demographics['PrepType'][index] == "MultiOme":
                if os.path.isfile(soloDoublet):
                    adata.obs["PrepType"] = ["multiOme"] * len(adata.obs_names)
                    csv_df = pd.read_csv(soloDoublet)
                    csv_index = list(csv_df.iloc[:, 0])
                    in_csv_dict = {index: index in csv_index for index in adata.obs_names}
                    adata.obs['SOLO'] = [str(in_csv_dict[name]) for name in adata.obs_names]
                    adata.obs['SOLO'] = adata.obs['SOLO'].replace({"True": 'doublet', 'False': 'singlet'})
                else:
                    print(demographics['SampleName'][index] + " does not have a SOLO barcode file.")
            else:
                print(demographics['SampleName'][index] + " is not a multiOme! No SOLO doublets to add.")
                adata.obs["SOLO"] = ['Not_measured'] * len(adata.obs_names)
                adata.obs["PrepType"] = ["snATAC"]  * len(adata.obs_names)
            print(type(adata.obs['SOLO'][1]))
            print(type(adata.obs['Amulet_doublet'][1]))
            adata.close()
            adata_input.close()
        else:
            print("No amulet file for " + demographics['SampleName'][index])

    else:
        print(demographics['SampleName'][index] + " does not have a post SOLO file!")


# In[5]:


doublet_col = demographics.columns.get_loc('snATAC_singlets')

for index, item in enumerate(demographics['fragment_files_snATAC']):
    output_file = '/home/levinsj/MultiOme/ATAC_data/Human/02_doublets_called/' + demographics['SampleName'][index]+ ".h5ad"
    if os.path.isfile(output_file):
        adata = sc.read_h5ad(output_file)
        print(adata)
        if 'Amulet_doublet' in adata.obs and 'singlet' in adata.obs['Amulet_doublet'].values:
            demographics.iloc[index, doublet_col] = adata.obs['Amulet_doublet'].value_counts()['singlet']
        else:
            print(demographics['SampleName'][index] + " does not doublets called correctly!")
    else:
        print(demographics['SampleName'][index] + " does not have a filtered h5ad file!")


# In[6]:


display(demographics)

demographics.to_csv('/home/levinsj/MultiOme/Sample_demographics.csv', index=False)

