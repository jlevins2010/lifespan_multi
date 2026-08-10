#!/usr/bin/env python
# coding: utf-8

# ### Time dependent Genes
# This script works to identify genes that change over time. It employs a pseudobulk approach and uses input from the cleaned integrated multiVI object. 

# ### Import Libraries

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata as ad
import scanpy as sc
import matplotlib as mpl
import os
import os.path
from IPython.display import display
import seaborn as sns

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


# ### scanpy variables

# In[2]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2
res = 2.0
n_features = 250000


# In[3]:


demographics = pd.read_csv('/home/levinsj/MultiOme/Sample_demographics.csv') 
demographics.index = demographics["SampleName"]


# ### SCANPY Settings

# In[4]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2


# In[5]:


inputSCVIfile = "/home/levinsj/MultiOme/MergedObjects/merged_postMultiVI_annotated_onlyGEX_human.h5ad"
output_dir = "/home/levinsj/MultiOme/Analysis/DEGs/human/"


# In[6]:


adata = sc.read_h5ad(inputSCVIfile, backed='r')


# #### Get cell types to test, so we can loop through them

# In[7]:


cell_types_to_test = ["PT","LOH","Podocyte","Endothelium","PEC","Stroma","Endothelium","DCT_CNT","IC"]


# ### Compare Fetal to Adult expression

# In[8]:


contrasts = ["type", "Adult", "Fetal"] ## increasing logCF means more expression over time
padj_threshold = 0.05
log2fc_threshold = 0.58 # Approximately 1.5-fold change in either direction (2^0.58 ~= 1.5)

for cell_type in cell_types_to_test:
    print(f"\n--- Processing Cell Type: {cell_type} ---")
    cell_type_subset = adata[adata.obs["cellType"] == cell_type].to_memory()
    sc.pp.filter_genes(cell_type_subset, min_counts=5)

    pseudo_bulk_samples = []
    unique_samples = cell_type_subset.obs["sample_final"].unique()

    if len(unique_samples) < 2:
        print(f"Warning: Not enough unique samples ({len(unique_samples)}) for DESeq2 analysis in {cell_type}. Skipping.")
        continue


    type_counts = cell_type_subset.obs.groupby('type')['sample_final'].nunique()
    if not (contrasts[1] in type_counts.index and contrasts[2] in type_counts.index and type_counts[contrasts[1]] >= 1 and type_counts[contrasts[2]] >= 1):
        print(f"Warning: Insufficient samples for {contrasts[1]} or {contrasts[2]} comparison in {cell_type}. Skipping.")
        continue

    for sample_id in unique_samples:
        sample_adata = cell_type_subset[cell_type_subset.obs["sample_final"] == sample_id]
        if hasattr(sample_adata.X, 'sum'): # Check if X is a sparse matrix that has a sum method
            summed_counts = sample_adata.X.sum(axis=0)
        else: # Assume it's a dense numpy array
            summed_counts = sample_adata.X.sum(axis=0)

        rep_adata = sc.AnnData(X=summed_counts.reshape(1, -1), var=cell_type_subset.var) # Use cell_type_subset.var to retain gene info
        rep_adata.obs_names = [sample_id]

        rep_adata.obs["type"] = sample_adata.obs["type"].iloc[0]
        pseudo_bulk_samples.append(rep_adata)

    pb_adata = sc.concat(pseudo_bulk_samples)
    print(pb_adata)
    counts_df = pd.DataFrame(pb_adata.X.astype(int), index=pb_adata.obs_names, columns=pb_adata.var_names)

    counts_df.columns = counts_df.columns.astype(str)

    if pb_adata.obs['type'].nunique() < 2:
        print(f"Warning: Only one unique 'type' found for {cell_type} after pseudo-bulking. Skipping DESeq2.")
        continue

    dds = DeseqDataSet(
        counts=counts_df,
        metadata=pb_adata.obs,
        design_factors="type",
    )

    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=contrasts) # Define contrast clearly
    stat_res.summary()
    ds_results = stat_res.results_df
    ds_results = ds_results.sort_values('stat', ascending=True)

    output_filepath = os.path.join(output_dir, f'{cell_type}_{contrasts[1]}_{contrasts[2]}_genes.csv')
    ds_results.to_csv(output_filepath, index=True)
    print(f"DESeq2 results saved to: {output_filepath}")


# ### Complare pediatric and fetal

# In[9]:


contrasts = ["type", "Pediatric", "Fetal"] #increasing logFC means more expression over time
padj_threshold = 0.05
log2fc_threshold = 0.58 # Approximately 1.5-fold change in either direction (2^0.58 ~= 1.5)

for cell_type in cell_types_to_test:
    print(f"\n--- Processing Cell Type: {cell_type} ---")
    cell_type_subset = adata[adata.obs["cellType"] == cell_type].to_memory()
    sc.pp.filter_genes(cell_type_subset, min_counts=5)

    pseudo_bulk_samples = []
    unique_samples = cell_type_subset.obs["sample_final"].unique()

    if len(unique_samples) < 2:
        print(f"Warning: Not enough unique samples ({len(unique_samples)}) for DESeq2 analysis in {cell_type}. Skipping.")
        continue


    type_counts = cell_type_subset.obs.groupby('type')['sample_final'].nunique()
    if not (contrasts[1] in type_counts.index and contrasts[2] in type_counts.index and type_counts[contrasts[1]] >= 1 and type_counts[contrasts[2]] >= 1):
        print(f"Warning: Insufficient samples for {contrasts[1]} or {contrasts[2]} comparison in {cell_type}. Skipping.")
        continue

    for sample_id in unique_samples:
        sample_adata = cell_type_subset[cell_type_subset.obs["sample_final"] == sample_id]
        if hasattr(sample_adata.X, 'sum'): # Check if X is a sparse matrix that has a sum method
            summed_counts = sample_adata.X.sum(axis=0)
        else: # Assume it's a dense numpy array
            summed_counts = sample_adata.X.sum(axis=0)

        rep_adata = sc.AnnData(X=summed_counts.reshape(1, -1), var=cell_type_subset.var) # Use cell_type_subset.var to retain gene info
        rep_adata.obs_names = [sample_id]

        rep_adata.obs["type"] = sample_adata.obs["type"].iloc[0]
        pseudo_bulk_samples.append(rep_adata)

    pb_adata = sc.concat(pseudo_bulk_samples)
    print(pb_adata)
    counts_df = pd.DataFrame(pb_adata.X.astype(int), index=pb_adata.obs_names, columns=pb_adata.var_names)

    counts_df.columns = counts_df.columns.astype(str)

    if pb_adata.obs['type'].nunique() < 2:
        print(f"Warning: Only one unique 'type' found for {cell_type} after pseudo-bulking. Skipping DESeq2.")
        continue

    dds = DeseqDataSet(
        counts=counts_df,
        metadata=pb_adata.obs,
        design_factors="type",
    )

    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=contrasts) # Define contrast clearly
    stat_res.summary()
    ds_results = stat_res.results_df
    ds_results = ds_results.sort_values('stat', ascending=True)

    output_filepath = os.path.join(output_dir, f'{cell_type}_{contrasts[1]}_{contrasts[2]}_genes.csv')
    ds_results.to_csv(output_filepath, index=True)
    print(f"DESeq2 results saved to: {output_filepath}")




# ### Compare Pediatric and Adult

# In[10]:


contrasts = ["type", "Adult", "Pediatric"] # increasing logFC means more expression over time
padj_threshold = 0.05
log2fc_threshold = 0.58 # Approximately 1.5-fold change in either direction (2^0.58 ~= 1.5)

for cell_type in cell_types_to_test:
    print(f"\n--- Processing Cell Type: {cell_type} ---")
    cell_type_subset = adata[adata.obs["cellType"] == cell_type].to_memory()
    sc.pp.filter_genes(cell_type_subset, min_counts=5)

    pseudo_bulk_samples = []
    unique_samples = cell_type_subset.obs["sample_final"].unique()

    if len(unique_samples) < 2:
        print(f"Warning: Not enough unique samples ({len(unique_samples)}) for DESeq2 analysis in {cell_type}. Skipping.")
        continue


    type_counts = cell_type_subset.obs.groupby('type')['sample_final'].nunique()
    if not (contrasts[1] in type_counts.index and contrasts[2] in type_counts.index and type_counts[contrasts[1]] >= 1 and type_counts[contrasts[2]] >= 1):
        print(f"Warning: Insufficient samples for {contrasts[1]} or {contrasts[2]} comparison in {cell_type}. Skipping.")
        continue

    for sample_id in unique_samples:
        sample_adata = cell_type_subset[cell_type_subset.obs["sample_final"] == sample_id]
        if hasattr(sample_adata.X, 'sum'): # Check if X is a sparse matrix that has a sum method
            summed_counts = sample_adata.X.sum(axis=0)
        else: # Assume it's a dense numpy array
            summed_counts = sample_adata.X.sum(axis=0)

        rep_adata = sc.AnnData(X=summed_counts.reshape(1, -1), var=cell_type_subset.var) # Use cell_type_subset.var to retain gene info
        rep_adata.obs_names = [sample_id]

        rep_adata.obs["type"] = sample_adata.obs["type"].iloc[0]
        pseudo_bulk_samples.append(rep_adata)

    pb_adata = sc.concat(pseudo_bulk_samples)
    print(pb_adata)
    counts_df = pd.DataFrame(pb_adata.X.astype(int), index=pb_adata.obs_names, columns=pb_adata.var_names)

    counts_df.columns = counts_df.columns.astype(str)

    if pb_adata.obs['type'].nunique() < 2:
        print(f"Warning: Only one unique 'type' found for {cell_type} after pseudo-bulking. Skipping DESeq2.")
        continue

    dds = DeseqDataSet(
        counts=counts_df,
        metadata=pb_adata.obs,
        design_factors="type",
    )

    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=contrasts) # Define contrast clearly
    stat_res.summary()
    ds_results = stat_res.results_df
    ds_results = ds_results.sort_values('stat', ascending=True)

    output_filepath = os.path.join(output_dir, f'{cell_type}_{contrasts[1]}_{contrasts[2]}_genes.csv')
    ds_results.to_csv(output_filepath, index=True)
    print(f"DESeq2 results saved to: {output_filepath}")



# ### Time dependent Genes

# In[11]:


cell_type = "NPC"
padj_threshold = 0.05

print(f"\n--- Processing Cell Type: {cell_type} ---")

mask = (~adata.obs['age'].isna()) & (adata.obs['age'] <= 1)
adata_cleaned = adata[mask, :].to_memory()
print(adata_cleaned.obs["sample_final"].value_counts())
cell_type_subset = adata_cleaned[adata_cleaned.obs["cellType"] == cell_type]
print(cell_type_subset.obs["sample_final"].value_counts())
sc.pp.filter_genes(cell_type_subset, min_counts=5)
print(cell_type_subset.obs["sample_final"].value_counts())

pseudo_bulk_samples = []
unique_samples = cell_type_subset.obs["sample_final"].unique()
type_counts = cell_type_subset.obs.groupby('type')['sample_final'].nunique()

for sample_id in unique_samples:
    sample_adata = cell_type_subset[cell_type_subset.obs["sample_final"] == sample_id]
    if hasattr(sample_adata.X, 'sum'): # Check if X is a sparse matrix that has a sum method
        summed_counts = sample_adata.X.sum(axis=0)
    else: # Assume it's a dense numpy array
        summed_counts = sample_adata.X.sum(axis=0)

    rep_adata = sc.AnnData(X=summed_counts.reshape(1, -1), var=cell_type_subset.var) # Use cell_type_subset.var to retain gene info
    rep_adata.obs_names = [sample_id]

    rep_adata.obs["age"] = sample_adata.obs["age"].iloc[0]
    rep_adata.obs["PrepType"] = sample_adata.obs["PrepType"].iloc[0]
    pseudo_bulk_samples.append(rep_adata)

pb_adata = sc.concat(pseudo_bulk_samples)
pb_adata.age = pb_adata.obs["age"].astype(float)


counts_df = pd.DataFrame(pb_adata.X.astype(int), index=pb_adata.obs_names, columns=pb_adata.var_names)
counts_df.columns = counts_df.columns.astype(str)

dds = DeseqDataSet(
    counts=counts_df,
    metadata=pb_adata.obs,
    design_factors='age',
    continuous_factors=['age']
)

dds.obs["age"] = dds.obs["age"].astype(float)
dds.deseq2()
stat_res = DeseqStats(dds, contrast=["age", "", ""])
stat_res.summary()
ds_results = stat_res.results_df
ds_results = ds_results.sort_values('padj', ascending = True)

output_filepath = os.path.join(output_dir, f'{cell_type}_timeDependent_genes.csv')
ds_results.to_csv(output_filepath, index=True)
print(f"DESeq2 results saved to: {output_filepath}")

