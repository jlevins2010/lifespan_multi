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
from glob import glob
from tqdm import tqdm
from IPython.display import display
import matplotlib as mpl


# ### scanpy and SNAPATAC2 settings

# In[2]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450

snap.__version__


# ### Import demographics

# In[3]:


demographics = pd.read_csv('/home/levinsj/MultiOme/Sample_demographics.csv') 
demographics.index = demographics["SampleName"]


# ### Create H5ad object

# In[4]:


snapATAC_col = demographics.columns.get_loc('snATAC_initialCells')

for index, item in enumerate(demographics['fragment_files_snATAC']):
    input_file = demographics['fragment_files_snATAC'][index]
    output_file = '/home/levinsj/MultiOme/ATAC_data/Human/raw_adata/' + demographics['SampleName'][index]+ ".h5ad"
    if os.path.isfile(input_file):
        if not os.path.isfile(output_file):
            adata = snap.pp.import_data(input_file, file=output_file, chrom_sizes=snap.genome.hg38, min_num_fragments=1000, sorted_by_barcode=False)
            snap.pp.add_tile_matrix(adata, bin_size=5000)
            snap.metrics.tsse(adata, snap.genome.hg38)
            demographics.iloc[index, snapATAC_col] = adata.n_obs
            adata.close()
            print(demographics['fragment_files_snATAC'][index] + " sample h5ad created!")
        else:
            print("h5ad already made for " + demographics['fragment_files_snATAC'][index])
    else:
        print(demographics['fragment_files_snATAC'][index] + " not found")


# In[5]:


display(demographics)

demographics.to_csv('/home/levinsj/MultiOme/Sample_demographics.csv', index=False)

