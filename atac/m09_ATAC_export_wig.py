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


input_anndataset = "/home/levinsj/MultiOme/ATAC_data/Human/03_integrated_w_doublets/postSCVI_final_data.h5ads/_dataset.h5ads"
filtered_barcodes = "/home/levinsj/MultiOme/Barcodes/snapATAC_Barcodes/type_celltype_anno_barcodes.csv"
output_anndataset = "/home/levinsj/MultiOme/ATAC_data/Human/04_SCVI/pairedSamples_final_data.h5ads/"
output_anndataset = "/home/levinsj/MultiOme/ATAC_data/Human/04_SCVI/postSCVI_ATAC_forBigWig.h5ad"


# ### load input file

# In[6]:


# Load the barcodes
snapBarcodes = pd.read_csv(filtered_barcodes, header=0) 
mapping_dict = {}

for _, row in snapBarcodes.iterrows():
    barcode_raw = str(row['barcode_orig']).strip()
    sample_id = str(row['sample']).strip()
    cell_type = str(row['annoCellType']).strip()

    # Match the exact format: HK2458_AAACCGGCATCATGTG-1
    mapping_dict[f"{sample_id}_{barcode_raw}"] = cell_type

print(f"Successfully built mapping dictionary with {len(mapping_dict)} unique cellular keys.")


# In[7]:


adataset = snap.read_dataset(filename=input_anndataset)
try:
    obs_series = pd.Series(adataset.obs_names)
    keepidx = obs_series.isin(mapping_dict.keys())
    adataset = adataset.subset(obs_indices=keepidx, out = output_anndataset)
    filtered_anndata = adataset[0]
    print("Dataset subsetted successfully!")
    print(filtered_anndata)
except Exception as e:
    print(f"An error occurred: {e}")


# In[8]:


try:
    mapped_list = [mapping_dict.get(str(cell).strip(), None) for cell in filtered_anndata.obs_names]
    filtered_anndata.obs["cellType_age"] = mapped_list
    valid_mappings = [x for x in mapped_list if x is not None]
    missing_count = mapped_list.count(None)

    print(f"Mapping completed successfully!")
    print(f"-> Total cells evaluated: {len(mapped_list)}")
    print(f"-> Successfully matched: {len(valid_mappings)}")
    print(f"-> Missing/Unmatched cells: {missing_count}")

    if len(valid_mappings) > 0:
        print("\nBreakdown of identified categories:")
        print(filtered_anndata.obs["cellType_age"].value_counts())

except Exception as e:
    print(f"An error occurred: {e}")


# In[9]:


try: 
    selections = ["Adult_PT","Pediatric_PT","Fetal_PT","Fetal_Int","Fetal_NPC","Adult_Stroma","Pediatric_Stroma", "Fetal_Stroma", "Adult_iPT","Fetal_Podocyte","Pediatric_Podocyte","Adult_Podocyte","Fetal_LOH","Adult_LOH","Pediatric_LOH", "Fetal_DCT_CNT","Adult_DCT_CNT","Pediatric_DCT_CNT"]
    snap.ex.export_coverage(filtered_anndata, groupby="cellType_age",selections = selections, suffix='.bw', normalization = "RPKM", include_for_norm = "/home/levinsj/Applications/hg38_tss_sites.bed",counting_strategy="insertion", bin_size=25, out_dir = "/home/levinsj/MultiOme/ATAC_data/Human/bigWig/")
    print("completed!!!")
except Exception as e:
    print(f"An error occurred: {e}")


# In[10]:


try:
    adataset[0].close()
    filtered_anndata.close()
except Exception as e:
    print(f"An error occurred: {e}")

