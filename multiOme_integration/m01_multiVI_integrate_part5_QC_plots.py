#!/usr/bin/env python
# coding: utf-8

# ## Further cleaning of the integrated object
# Cleaning and annotation of the multiVI object. It also subsets to just paired cells. Run on CUDA.

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
from pathlib import Path


# In[2]:


colors = {"DCT_CNT": "black",
               "Endothelium": "#7ae031",
               "Podocyte": "#b56a07", 
               "Stroma": "#794b82",
               "NPC": "#ff8000", 
               "PT": "#ff00d4", 
               "Int": "#698cff",
               "IC": "#191985", 
               "PEC": "#ff0011", 
               "LOH": "#235e00",
               "Immune Cells": '#757575',
               "iPT":"#f78e52"
         }



# In[3]:


sc.set_figure_params(
    dpi=80,  
    dpi_save=300,  
    frameon=False, 
    vector_friendly=True, 
    format="pdf",  
)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['pdf.fonttype'] = 42 
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.size'] = 7
plt.rcParams['axes.titlesize'] = 8
plt.rcParams['axes.labelsize'] = 7
plt.rcParams['xtick.labelsize'] = 6
plt.rcParams['ytick.labelsize'] = 6

sc.settings.figdir = "/home/levinsj/MultiOme/Analysis/pythonNotebooks/pythonPlots"
fig_savedir = Path("/home/levinsj/MultiOme/Analysis/pythonNotebooks/pythonPlots")


# ### scanpy variables

# In[4]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2
res = 2.0
n_features = 250000


# ### SCVI import

# In[5]:


batch_size = 32
max_epochs = 2000
earlyStopping = True
learning_rate = 0.0005
scvi_layers = 3
scvi_latent = 30
fraction_expressed = 0.01

multiVI_LATENT_KEY = "X_multivi"
multiVI_CLUSTERS_KEY = "clusters_multivi"


# ### SCANPY Settings

# In[6]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2


# ### Helpful function

# In[7]:


def keep_matching_elements(my_list, adata_var_names):
  return [element for element in my_list if element in adata_var_names]


# In[8]:


annotatedFile = "/home/levinsj/MultiOme/MergedObjects/merged_postMultiVI_annotated_final_human.h5ad"


# In[9]:


adata = sc.read_h5ad(annotatedFile)

#sc.pl.umap(adata, color = multiVI_CLUSTERS_KEY, legend_loc='on data', frameon = False, legend_fontsize=10, legend_fontoutline=2)
sc.pl.umap(adata, color = "PrepType", save="_prepType_annotated.pdf", show=True)
sc.pl.umap(adata, color = "study", save="_study_annotated.pdf", show=True)
sc.pl.umap(adata, color = "type", save="_type_annotated.pdf", show=True)
sc.pl.umap(adata, color = "sample_final", save="_sampleFinal_annotated.pdf", show=True)


# In[10]:


sc.pl.umap(adata, color = "n_genes_by_counts", save="_n_genes_by_counts_annotated.pdf", show=True)
sc.pl.umap(adata, color = "total_counts", save="_total_counts_annotated.pdf", show=True)
sc.pl.umap(adata, color = "pct_counts_mt", save="_pct_counts_mt.pdf", show=True)


# In[11]:


for i in adata.obs["type"].unique():

    var_names = ['UNCX', 'ITGA8', 'LHX1', 'JAG1', 'PTPRO', 'MAFB','WT1',"CFH","VCAM1","HAVCR1", 'CUBN', 'SLC13A1', 'SLC12A1', 'UMOD', 'SLC12A3', 'GATA3', 'SLC26A4', 'ATP6V1B1', 'PLVAP', 'EMCN', 'COL1A1', 'COL3A1', 'DOCK2', 'PTPRC']
    order = ["NPC","Int","Podocyte","PEC","iPT", "PT","LOH","DCT_CNT","IC","Endothelium","Stroma","Immune Cells"]
    print(i)
    sc.pl.dotplot(
            adata[adata.obs["type"] == i],
            var_names=keep_matching_elements(var_names, adata.var_names),
            groupby='cellType',
            categories_order = order,
            standard_scale='var',  # Scale gene expression values per gene to [0, 1]
            cmap='Blues',        # Colormap for expression intensity
            dot_max=0.6,           # Max relative dot size
            dot_min=0.0,           # Min relative dot size
            figsize=(8, 6), 
            dendrogram=False,
            log=True, save=f"_{i}_dotplot.pdf", show=True)


# In[12]:


df = adata.obs[multiVI_CLUSTERS_KEY].value_counts().to_frame()
df.index = adata.obs[multiVI_CLUSTERS_KEY].value_counts().index

df['sampleBreakdown'] = object
df['PrepType'] = object
df['modality'] = object

samples = list(set(adata.obs["sample_final"].values.ravel().tolist()))

for i in df.index:
    counts_sample = []
    prep_sample = []
    modality_sample = []
    for j in samples:
        counts_sample.append(adata.obs[(adata.obs['sample_final'] == j) & (adata.obs[multiVI_CLUSTERS_KEY] == i)].shape[0]/adata.obs[(adata.obs[multiVI_CLUSTERS_KEY] == i)].shape[0])
    prep_sample.append(adata.obs[(adata.obs['PrepType'] == "multiOme") & (adata.obs[multiVI_CLUSTERS_KEY] == i)].shape[0]/adata.obs[(adata.obs[multiVI_CLUSTERS_KEY] == i)].shape[0])
    modality_sample.append(adata.obs[(adata.obs['modality'] == "paired") & (adata.obs[multiVI_CLUSTERS_KEY] == i)].shape[0]/adata.obs[(adata.obs[multiVI_CLUSTERS_KEY] == i)].shape[0])
    df.at[i,'sampleBreakdown'] = counts_sample
    df.at[i,'PrepType'] = prep_sample
    df.at[i,'modality'] = modality_sample


# In[13]:


width = 0.9

# Loop over each cell type and create a stacked bar plot of the "sample" category
for i, cell_type in enumerate(df.index):
    ind = np.arange(1)
    S1 = df.sampleBreakdown[i]

    # generate a palette with all samples (scramble colors)
    num_colors = len(samples)
    cmap = plt.cm.get_cmap('hsv', num_colors)
    palette_indices = np.random.permutation(num_colors)
    palette = [cmap(i) for i in palette_indices]

    bottom = 0
    for j in range(len(S1)):
          plt.bar(i, S1[j],bottom = bottom, width = width, color = palette[j])
          bottom = bottom + S1[j]
    # Show the plot
plt.tight_layout()

plt.show()


# In[14]:


df = adata.obs["cellType"].value_counts().to_frame()
df.index = adata.obs["cellType"].value_counts().index

df['sampleBreakdown'] = object
df['sample'] = object
df['composition'] = object

samples = list(["Adult","Fetal","Pediatric"])

for i in df.index:
    counts_sample = []
    for j in samples:
        counts_sample.append(adata.obs[(adata.obs['type'] == j) & (adata.obs['cellType'] == i)].shape[0])
    df.at[i,'sampleBreakdown'] = counts_sample

print(df)    
cell_types = df.index
fig, ax = plt.subplots(figsize=(12, 6))
width = 0.9

tab10_colors_mpl = plt.get_cmap('tab10').colors[:3]
colors = [
    tab10_colors_mpl[0],  # Blue for Adult
    tab10_colors_mpl[1],  # Orange for Fetal
    tab10_colors_mpl[2]   # Green for Pediatric
]
print(colors)

# Loop through each sample and plot the bar segments
bottoms = [0] * len(cell_types)
for j, sample in enumerate(samples):
    # Get the counts for the current sample across all cell types
    sample_counts = [df.sampleBreakdown[i][j] for i in range(len(cell_types))]

    # Plot the bar segment for the current sample
    ax.bar(np.arange(len(cell_types)), sample_counts,
           bottom=bottoms,
           label=sample,
           color=colors[j],
           width=width)

    # Update the bottom for the next sample
    bottoms = [bottoms[k] + sample_counts[k] for k in range(len(cell_types))]

# Add labels and legend
ax.set_xticks(np.arange(len(cell_types)))
ax.set_xticklabels(cell_types, rotation=90) # Rotate for long labels
ax.set_ylabel('Cell Count')
ax.set_title('Cell Type Composition by Type')

plt.tight_layout()

figure_name = "CellType_composition_QCplot.pdf"
save_path = fig_savedir / figure_name
plt.savefig(save_path, format='pdf', bbox_inches='tight')
plt.show()


# In[15]:


df = adata.obs["cellType"].value_counts().to_frame()
df.index = adata.obs["cellType"].value_counts().index

df['sampleBreakdown'] = object
df['sample'] = object
df['composition'] = object

samples = list(["multiOme","snATAC","snRNA"])

for i in df.index:
    counts_sample = []
    for j in samples:
        counts_sample.append(adata.obs[(adata.obs['PrepType'] == j) & (adata.obs['cellType'] == i)].shape[0])
    df.at[i,'sampleBreakdown'] = counts_sample

print(df)    
cell_types = df.index
fig, ax = plt.subplots(figsize=(12, 6))
width = 0.9

tab10_colors_mpl = plt.get_cmap('tab10').colors[:3]
colors = [
    tab10_colors_mpl[0],  # Blue for Adult
    tab10_colors_mpl[1],  # Orange for Fetal
    tab10_colors_mpl[2]   # Green for Pediatric
]
print(colors)

# Loop through each sample and plot the bar segments
bottoms = [0] * len(cell_types)
for j, sample in enumerate(samples):
    # Get the counts for the current sample across all cell types
    sample_counts = [df.sampleBreakdown[i][j] for i in range(len(cell_types))]

    # Plot the bar segment for the current sample
    ax.bar(np.arange(len(cell_types)), sample_counts,
           bottom=bottoms,
           label=sample,
           color=colors[j],
           width=width)

    # Update the bottom for the next sample
    bottoms = [bottoms[k] + sample_counts[k] for k in range(len(cell_types))]

# Add labels and legend
ax.set_xticks(np.arange(len(cell_types)))
ax.set_xticklabels(cell_types, rotation=90) # Rotate for long labels
ax.set_ylabel('Cell Count')
ax.set_title('Cell Type Composition by Type')

plt.tight_layout()
figure_name = "PrepType_composition_QCplot.pdf"
save_path = fig_savedir / figure_name
plt.savefig(save_path, format='pdf', bbox_inches='tight')
plt.show()


# In[16]:


df = adata.obs["cellType"].value_counts().to_frame()
df.index = adata.obs["cellType"].value_counts().index

df['sampleBreakdown'] = object
df['sample'] = object
df['composition'] = object

samples = list(set(adata.obs["sample_final"].values.ravel().tolist()))

for i in df.index:
    counts_sample = []
    for j in samples:
        counts_sample.append(adata.obs[(adata.obs['sample_final'] == j) & (adata.obs['cellType'] == i)].shape[0])
    df.at[i,'sampleBreakdown'] = counts_sample

print(df)    

width = 0.9
fig, ax = plt.subplots(figsize=(12, 6))
width = 0.9

# Loop over each cell type and create a stacked bar plot of the "sample" category
for i, cell_type in enumerate(df.index):
    ind = np.arange(1)
    S1 = df.sampleBreakdown[i]

    # generate a palette with all samples (scramble colors)
    num_colors = len(samples)
    cmap = plt.cm.get_cmap('hsv', num_colors)
    palette_indices = np.random.permutation(num_colors)
    palette = [cmap(i) for i in palette_indices]

    bottom = 0
    for j in range(len(S1)):
          ax.bar(i, S1[j],bottom = bottom, width = width, color = palette[j])
          bottom = bottom + S1[j]
    # Show the plot
ax.set_xticks(np.arange(len(cell_types)))
ax.set_xticklabels(cell_types, rotation=90) # Rotate for long labels
ax.set_ylabel('Cell Count')
ax.set_title('Cell Type Composition by Sample')

plt.tight_layout()
figure_name = "sample_composition_QCplot.pdf"
save_path = fig_savedir / figure_name
plt.savefig(save_path, format='pdf', bbox_inches='tight')
plt.show()

