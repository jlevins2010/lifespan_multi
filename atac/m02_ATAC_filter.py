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
import plotly.graph_objects as go
import shutil


# ### Scanpy and SNAPATAC2 settings

# In[2]:


snap.__version__

sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450

black_list_file = '/home/levinsj/Applications/hg38-blacklist.bed'


# ### Demographics

# In[3]:


demographics = pd.read_csv('/home/levinsj/MultiOme/Sample_demographics.csv') 
demographics.index = demographics["SampleName"]


# ### filter the adata object for QC and export barcodes of all cells called

# In[4]:


snapATAC_col = demographics.columns.get_loc('snATAC_initialCells')
postFilter_col = demographics.columns.get_loc('snATAC_post_filter')

for index, item in enumerate(demographics['fragment_files_snATAC']):
    input_file = '/home/levinsj/MultiOme/ATAC_data/Human/raw_adata/' + demographics['SampleName'][index]+ ".h5ad"
    output_file = '/home/levinsj/MultiOme/ATAC_data/Human/01_filtered/' + demographics['SampleName'][index]+ ".h5ad"
    barcodes_output_file = "/home/levinsj/MultiOme/Barcodes/snapATAC_Barcodes/" +demographics['SampleName'][index]+"_snapBarcodes.csv"
    ### import the sample level cutoffs
    min_frags = demographics['min_frags'][index]
    min_tsse = demographics['min_tsse'][index]
    max_frags = demographics['max_frags'][index]
    if os.path.isfile(input_file):
        if not os.path.isfile(output_file):
            try:
                adata_input = snap.read(filename=input_file)
                print(demographics['SampleName'][index])
                fig_pre_filter = snap.pl.tsse(adata_input, interactive=False) # show=False prevents display
                display(fig_pre_filter)
                adata = adata_input.copy(output_file)
                demographics.iloc[index, snapATAC_col] = adata.n_obs          
                # filter ATAC samples
                snap.pp.filter_cells(adata, min_counts=min_frags, min_tsse=min_tsse, max_counts=max_frags)
                snap.pp.select_features(adata, blacklist=black_list_file)
                fig_post_filter = snap.pl.tsse(adata, interactive=False)
                display(fig_post_filter)
                demographics.iloc[index, postFilter_col] = adata.n_obs
                print(demographics['SampleName'][index] + " sample h5ad created!")
                # export all called indecies
                index_df = pd.DataFrame(adata.obs_names)
                index_df.to_csv(barcodes_output_file, index=False, header=False)  
                adata.close()
                adata_input.close()
            except Exception as e:
                print("An unexpected error occurred:", e)
        else:
            print("h5ad already made for " + demographics['SampleName'][index])
    else:
        print(demographics['SampleName'][index] + " not found")


# ### Create a singlecell.csv file for each sample for amulet
# need to move the per-sample barcodes from the cellranger output to correct file and rename as specified below.

# In[5]:


for index, item in enumerate(demographics['fragment_files_snATAC']):
    try:
        input_file = demographics['fragment_files_snATAC'][index]
        parent_dir = os.path.dirname(input_file)
        perSampleBarcodePath = os.path.join(parent_dir, "per_barcode_metrics.csv")
        destination_path = '/home/levinsj/MultiOme/Barcodes/perSampleBarcodes/' + demographics['SampleName'][index]+ "_per_sample_barcodes.csv"
        shutil.copy2(perSampleBarcodePath, destination_path)
        print(destination_path + "created")
    except FileNotFoundError:
        print(f"Error: Source file not found.")


# In[6]:


for index, item in enumerate(demographics['fragment_files_snATAC']):
    per_barcode_file = '/home/levinsj/MultiOme/Barcodes/perSampleBarcodes/' + demographics['SampleName'][index]+ "_per_sample_barcodes.csv"
    barcodes_output_file = "/home/levinsj/MultiOme/Barcodes/snapATAC_Barcodes/" +demographics['SampleName'][index]+"_snapBarcodes.csv"
    amulet_imput = "/home/levinsj/MultiOme/Barcodes/amulet_input/" +demographics['SampleName'][index]+"_singlecell.csv"
    if os.path.isfile(barcodes_output_file):
        if os.path.isfile(per_barcode_file):
            if not os.path.isfile(amulet_imput):
                # get all barcodes called as cells from SnapATAC2
                snapBarcodes = pd.read_csv(barcodes_output_file) 
                calledBarcodes = snapBarcodes.iloc[:, 0].values
                # check if these barcodes are in the the persampleBarcodes file
                perSampleBarcodes = pd.read_csv(per_barcode_file) 
                perSampleBarcodes["is__cell_barcode"] = perSampleBarcodes.apply(lambda row: 1 if row["barcode"] in calledBarcodes else 0, axis=1)
                perSampleBarcodes = perSampleBarcodes[["barcode","is__cell_barcode"]]
                perSampleBarcodes.to_csv(amulet_imput, index=False)
                print(demographics['SampleName'][index]+" is ready for Amulet.")
            else: 
                print(demographics['SampleName'][index]+" was already run.")
        else:
            print(demographics['SampleName'][index]+" does not have the correct per_barcode_file from cellranger or is already written.")
    else:
        print(demographics['SampleName'][index]+" does not have a SnapAtac2 barcodes files!")


# In[7]:


display(demographics)
demographics.to_csv('/home/levinsj/MultiOme/Sample_demographics.csv', index=False)

