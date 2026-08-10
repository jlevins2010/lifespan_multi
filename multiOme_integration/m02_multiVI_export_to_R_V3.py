#!/usr/bin/env python
# coding: utf-8

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
from scipy import sparse
from scipy.sparse import csr_matrix
import scipy.io as sio


# ### scanpy variables

# In[2]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2
res = 2.0
n_features = 250000


# ### SCANPY Settings

# In[3]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2
np.random.seed(4)


# ### demographics

# In[4]:


outputSCVIfile_paired_only = "/home/levinsj/MultiOme/MergedObjects/merged_postMultiVI_annotated_onlyPaired_human.h5ad"

adata_paired_GEX_path_forGWAS = "/home/levinsj/MultiOme/MergedObjects/merged_postMultiVI_annotated_onlyPaired_human_GEX_GWAS.h5ad"
adata_paired_ATAC_path_forGWAS = "/home/levinsj/MultiOme/MergedObjects/merged_postMultiVI_annotated_onlyPaired_human_ATAC_GWAS.h5ad"

adata_paired_GEX_path_PT = "/home/levinsj/MultiOme/MergedObjects/merged_postMultiVI_annotated_onlyPaired_human_GEX_forNPC.h5ad"
adata_paired_ATAC_path_PT = "/home/levinsj/MultiOme/MergedObjects/merged_postMultiVI_annotated_onlyPaired_human_ATAC_forNPC.h5ad"

all_outputs = [adata_paired_GEX_path_forGWAS, adata_paired_ATAC_path_forGWAS, adata_paired_GEX_path_PT, adata_paired_ATAC_path_PT]
all_names = ["GEX_GWAS", "ATAC_GWAS", "GEX_PT", "ATAC_PT"]


# In[5]:


obs_to_keep = ["sample_final", "cellType", "type","combined_type"]


# In[6]:


### Functions


# In[7]:


def downsample_adata(adata_obj, group_col, n_samples):
    sampled_indices = (
        adata_obj.obs.groupby(group_col)
        .apply(lambda x: x.sample(n=n_samples, random_state=42))
        .index.get_level_values(1)
    )
    return adata_obj[sampled_indices, :].copy()


def process_and_export_modality(adata_obj, modality_id, output_path, cols_to_keep):
    adata_mod = adata_obj[:, adata_obj.var.modality == modality_id].copy()
    adata_mod.obs_names = (
        adata_mod.obs["sample_final"].astype(str) + "_" + adata_mod.obs_names
    )

    adata_mod.obs = adata_mod.obs[cols_to_keep].astype(str)
    for col in adata_mod.var.columns:
        if adata_mod.var[col].dtype.name == "category":
            adata_mod.var[col] = adata_mod.var[col].astype(str)

    adata_mod.obsp = {}
    adata_mod.obsm = {}
    adata_mod.layers = {}
    adata_mod.uns = {}

    if not isinstance(adata_mod.X, csr_matrix):
        adata_mod.X = csr_matrix(adata_mod.X)

    adata_mod.write_h5ad(filename=output_path)
    return adata_mod


# In[8]:


adata_paired = sc.read_h5ad(outputSCVIfile_paired_only)
adata_paired.obs_names_make_unique()

adata_paired.obs.loc[adata_paired.obs["sample_final"] == "HK3558", "type"] = (
    "Fetal"
)

adata_paired.obs["combined_type"] = (
    adata_paired.obs["type"].astype(str)
    + "_"
    + adata_paired.obs["cellType"].astype(str)
)

type_counts = adata_paired.obs["type"].value_counts()
min_count_gwas = type_counts.min()
print(f"Downsampling GWAS groups to minimum size: {min_count_gwas}")

adata_gwas_downsampled = downsample_adata(
    adata_paired, "type", n_samples=min_count_gwas
)

adata_gwas_GEX = process_and_export_modality(
    adata_gwas_downsampled, "0", adata_paired_GEX_path_forGWAS, obs_to_keep
)
adata_gwas_ATAC = process_and_export_modality(
    adata_gwas_downsampled, "1", adata_paired_ATAC_path_forGWAS, obs_to_keep
)


# In[9]:


adata_paired = sc.read_h5ad(outputSCVIfile_paired_only)
adata_paired.obs_names_make_unique()

# Concatenate columns into a single identifier
adata_paired.obs["combined_type"] = (
    adata_paired.obs["type"].astype(str)
    + "_"
    + adata_paired.obs["cellType"].astype(str)
)

# Target specific subsets
target_groups = ["Adult_PT", "Pediatric_PT", "Fetal_PT", "Fetal_Int", "Fetal_NPC"]
adata_pt_subset = adata_paired[
    adata_paired.obs["combined_type"].isin(target_groups)
].copy()

# Determine minimum group size within the target categories
min_count_pt = adata_pt_subset.obs["combined_type"].value_counts().min()
print(f"Downsampling PT/NPC groups to minimum size: {min_count_pt}")

# Downsample via our function
adata_pt_downsampled = downsample_adata(
    adata_pt_subset, "combined_type", n_samples=min_count_pt
)

# Print verification of perfect equal sizing
print("Final Cell Counts per Group:")
print(adata_pt_downsampled.obs["combined_type"].value_counts())

# Process and save GEX and ATAC for PT/NPC
adata_pt_GEX = process_and_export_modality(
    adata_pt_downsampled, "0", adata_paired_GEX_path_PT, obs_to_keep
)
adata_pt_ATAC = process_and_export_modality(
    adata_pt_downsampled, "1", adata_paired_ATAC_path_PT, obs_to_keep
)
print("PT/NPC Export Complete.")


# In[10]:


output_dir = "/home/levinsj/MultiOme/MergedObjects/R_objects/conversion_files/"
os.makedirs(output_dir, exist_ok=True)

for j, i in enumerate(all_outputs):
    print(f"\n--- Processing Export for Modality Group: {all_names[j]} ---")
    print(f"Loading AnnData object from: {i}")
    adata = sc.read_h5ad(i)

    # 1. Matrix Market Coordinate format (.T for GenexCell R format)
    matrix_path = os.path.join(output_dir, f"{all_names[j]}_matrix.mtx")
    print(f"Writing matrix chunk to: {matrix_path}")
    sio.mmwrite(matrix_path, adata.X.T)

    # 2. Unique Cell Barcodes
    barcodes_path = os.path.join(output_dir, f"{all_names[j]}_barcodes.tsv")
    print(f"Writing barcodes to: {barcodes_path}")
    with open(barcodes_path, "w") as f:
        f.write("\n".join(adata.obs_names.tolist()) + "\n")

    # 3. Unique Gene/Peak Features
    genes_path = os.path.join(output_dir, f"{all_names[j]}_features.tsv")
    print(f"Writing features to: {genes_path}")
    with open(genes_path, "w") as f:
        f.write("\n".join(adata.var_names.tolist()) + "\n")

    # 4. Clean Metadata DataFrame
    print("Exporting clean cell metadata...")
    metadata_df = adata.obs.copy()
    for col in metadata_df.columns:
        if (
            metadata_df[col].dtype.name == "category"
            or metadata_df[col].dtype.name == "object"
        ):
            metadata_df[col] = metadata_df[col].astype(str).fillna("Unknown")

    metadata_out_path = os.path.join(output_dir, f"{all_names[j]}_metadata.csv")
    metadata_df.to_csv(metadata_out_path)

print(
    f"\nSuccess! All pristine primitive files have been cleanly isolated to:\n{output_dir}"
)

