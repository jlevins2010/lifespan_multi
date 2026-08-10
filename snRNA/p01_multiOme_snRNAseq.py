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
import seaborn as sns
from scipy.stats import median_abs_deviation
import os.path
from IPython.display import display


# ### Scanpy Settings

# In[2]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450

minGenes = 200 # use only cells with at least 200 genes
minCells = 3 # use only genes expressed in at last 3 cells
mtThresh = 50


# ### Demographics

# In[3]:


demographics = pd.read_csv('/home/levinsj/MultiOme/Sample_demographics.csv') 
demographics.index = demographics["SampleName"]
display(demographics)


# ### Create h5ad for each sample

# In[4]:


for index, item in enumerate(demographics['sampleRawH5Matricies_snRNA']):
    snRNA_matrix_file = demographics['sampleRawH5Matricies_snRNA'][index]
    output_file = "/home/levinsj/MultiOme/CellBenderCorrected_snRNA/01_rawCellBender/h5ad/"+demographics['SampleName'][index]+"_cellBender.h5ad"
    if os.path.isfile(snRNA_matrix_file):
        if not os.path.isfile(output_file):
            adata = sc.read_10x_h5(snRNA_matrix_file, genome="GRCh38")
            adata.write_h5ad(output_file)
        else:
            print("h5ad already made for " + demographics['SampleName'][index])
    else:
        print(demographics['SampleName'][index] + " h5ad not found")


# ### Preprocess samples

# In[5]:


def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


# In[6]:


cellBenderCell_col = demographics.columns.get_loc('snRNA_cellBender_cells')
post_filter_col = demographics.columns.get_loc('snRNA_postInitialFilter_cells')

for index, item in enumerate(demographics['sampleRawH5Matricies_snRNA']):
    input_file = "/home/levinsj/MultiOme/CellBenderCorrected_snRNA/01_rawCellBender/h5ad/"+demographics['SampleName'][index]+"_cellBender.h5ad"
    output_file = "/home/levinsj/MultiOme/CellBenderCorrected_snRNA/02_initialQC/"+demographics['SampleName'][index]+"_initialQC.h5ad"
    if os.path.isfile(input_file):
        if not os.path.isfile(output_file):
            adata = sc.read_h5ad(input_file)
            demographics.iloc[index, cellBenderCell_col] = adata.n_obs

            adata.var_names_make_unique()
            adata.var["mt"] = adata.var_names.str.startswith("MT-")
            adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
            adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

            adata.obs["sample"] = demographics['SampleName'][index]
            adata.obs["study"] = demographics['Site'][index]
            adata.obs["type"] = demographics['Type'][index]
            adata.layers["counts"] = adata.X.copy()
            sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
            adata.obs["outlier"] = (
            is_outlier(adata, "log1p_total_counts", 5)
            | is_outlier(adata, "log1p_n_genes_by_counts", 5)
            | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
            )

            adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
            adata.obs["pct_counts_mt"] > mtThresh)

            adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

            adata.layers["counts"] = adata.X.copy()
            sc.pp.normalize_per_cell(adata)
            sc.pp.log1p(adata)
            #sc.pp.filter_genes(adata, min_cells=minCells) ## removed 7/15/25
            sc.pp.filter_cells(adata, min_genes=minGenes)
            adata.write(output_file)
            demographics.iloc[index, post_filter_col] = adata.n_obs
        else:
            print(demographics['SampleName'][index] + " is aready pre-processed!")
    else:
        print(demographics['SampleName'][index] + " does not have a cellBender h5ad file!")


# In[7]:


display(demographics)

demographics.to_csv('/home/levinsj/MultiOme/Sample_demographics.csv', index=False)

