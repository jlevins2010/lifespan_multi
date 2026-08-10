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
from matplotlib import cm # Import colormap module

from pygam import GAM, s, te, PoissonGAM, GammaGAM, LinearGAM, l, LogisticGAM # s() for smoothing spline, l() for linear term
from pygam.utils import OptimizationError

import seaborn as sns
import scipy


# In[2]:


colors = {"DCT_CNT": "black",
               "Endothelium": "#7ae031",
               "Podocyte": "#ad9c00", 
               "Stroma": "#794b82",
               "NPC": "#ff8000", 
               "PT": "#ff00d4", 
               "Int": "#698cff",
               "IC": "#191985", 
               "PEC": "#ff0011", 
               "LOH": "#235e00",
               "Immune Cells": '#757575',
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


# ### scanpy variables

# In[4]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2
res = 2.0
n_features = 250000


# In[5]:


demographics = pd.read_csv('/home/levinsj/MultiOme/Sample_demographics.csv') 
demographics.index = demographics["SampleName"]

keys = demographics['SampleName']
age = demographics['Age']
gAge = demographics['gAge']

# Use zip to pair them up and then convert to a dictionary
age_dict = dict(zip(keys, age))
gAge_dict = dict(zip(keys, gAge))

age_dict['HK2801'] = '45'
age_dict['HK3032'] = '61'
age_dict['HK2826'] = '68'
age_dict['HK2804'] = '49'
age_dict['HK3101'] = '48'
age_dict['HK2823'] = '20'
gAge_dict['HK2801'] = 'Not documented'
gAge_dict['HK3032'] = 'Not documented'
gAge_dict['HK2826'] = 'Not documented'
gAge_dict['HK2804'] = 'Not documented'
gAge_dict['HK3101'] = 'Not documented'
gAge_dict['HK2823'] = 'Not documented'


# ### SCANPY Settings

# In[6]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2


# In[7]:


inputSCVIfile = "/home/levinsj/MultiOme/MergedObjects/merged_postMultiVI_annotated_onlyGEX_human.h5ad"
output_dir = "/home/levinsj/MultiOme/Analysis/DEGs/"


# In[8]:


adata = sc.read_h5ad(inputSCVIfile)
adata.X = adata.layers["counts"]

# 2. Map metadata
adata.obs["age_mapped"] = adata.obs["sample_final"].map(age_dict)
adata.obs["gestational_age_mapped"] = adata.obs["sample_final"].map(gAge_dict)

non_nan_age_mask = adata.obs["age_mapped"].notna()
adata_cleaned = adata[non_nan_age_mask, :].copy()

# 4. Safely convert mapped age to numeric
adata_cleaned.obs["age"] = pd.to_numeric(
    adata_cleaned.obs["age_mapped"], errors="coerce"
)

ga_series = adata_cleaned.obs["gestational_age_mapped"].astype(str)
ga_numeric = pd.to_numeric(ga_series.str.extract(r"(\d+)", expand=False), errors="coerce")

# 6. Apply gestational age formula
adata_cleaned.obs["gestational_age"] = 0 - ((39 - ga_numeric) / 52)

# 7. Update age where gestational age exists
has_gestational_age = adata_cleaned.obs["gestational_age"].notna()
adata_cleaned.obs.loc[has_gestational_age, "age"] = adata_cleaned.obs.loc[
    has_gestational_age, "gestational_age"
]

# 8. Clean up intermediate temporary columns (Optional but recommended)
adata_cleaned.obs.drop(
    columns=["age_mapped", "gestational_age_mapped"], inplace=True
)

# 9. Save cleaned AnnData
adata_cleaned.write_h5ad(
    "/home/levinsj/MultiOme/MergedObjects/merged_postMultiVI_annotated_onlyGEX_human_cleaned.h5ad"
)


# In[9]:


age_data = adata_cleaned.obs['age'].dropna()

print(f"Number of age values after cleaning NaNs: {len(age_data)}")

plt.figure(figsize=(10, 6)) 
sns.histplot(age_data, bins='auto', color='skyblue', edgecolor='black')

plt.title('Distribution of Age', fontsize=16, fontweight='bold')
plt.xlabel('Age', fontsize=12)
plt.ylabel('Number of Cells', fontsize=12)

plt.grid(axis='y', alpha=0.75, linestyle='--')

plt.tight_layout() 
plt.show()


# In[10]:


sc.pl.violin(adata_cleaned, keys = ['UMOD'],  size = 1, groupby = 'type', rotation= 45, layer = "counts", use_raw = False)


# #### Get cell types to test, so we can loop through them

# In[11]:


SAMPLE_COL = "sample_final"  
CELL_TYPE_COL = "cellType"
AGE_COL = "age"

TARGET_CELL_TYPES = ["PT", "LOH", "DCT_CNT", "Stroma", "Endothelium", "Podocyte", "PEC", "Immune Cells", "iPT"] 

# 1. Base sample summary from cleaned dataset
sample_summary = adata_cleaned.obs.groupby(SAMPLE_COL).agg(
    total_cells=(CELL_TYPE_COL, 'size'),
    age=(AGE_COL, 'mean')
).reset_index()

# 2. Add counts & fractions per cell type using adata_cleaned (NOT adata)
for cell_type in TARGET_CELL_TYPES:
    indicator_col = f'is_{cell_type}'
    adata_cleaned.obs[indicator_col] = (adata_cleaned.obs[CELL_TYPE_COL] == cell_type).astype(int)

    # Calculate sum per sample
    cell_counts = adata_cleaned.obs.groupby(SAMPLE_COL)[indicator_col].sum().reset_index()
    cell_counts.rename(columns={indicator_col: f'count_{cell_type}'}, inplace=True)

    # Merge into summary table
    sample_summary = pd.merge(sample_summary, cell_counts, on=SAMPLE_COL)

    # Compute proportion
    sample_summary[f'fraction_{cell_type}'] = sample_summary[f'count_{cell_type}'] / sample_summary['total_cells']

print("Sample Summary DataFrame with Fractions:")
print(sample_summary.head())


# In[12]:


N_SPLINES = 6
LAMBDA = 2


# In[13]:


sns.set_style("whitegrid")

fig, axes = plt.subplots(
    nrows=1, 
    ncols=len(TARGET_CELL_TYPES), 
    figsize=(2 * len(TARGET_CELL_TYPES), 2.5) # Scaled down figure size to fit 6-8pt fonts nicely
)

if len(TARGET_CELL_TYPES) == 1:
    axes = [axes]
else:
    axes = axes.flatten()

X_age = sample_summary[[AGE_COL]].values 
age_range_pred = np.linspace(X_age.min(), X_age.max(), 100).reshape(-1, 1)

def inverse_logit(x):
    return 1 / (1 + np.exp(-x))

EPSILON = 1e-6 

for i, cell_type in enumerate(TARGET_CELL_TYPES):
    ax = axes[i]
    line_color = colors.get(cell_type, 'gray')

    # 1. Transformed response
    y_fraction = sample_summary[f'fraction_{cell_type}'].values
    y_clipped = np.clip(y_fraction, EPSILON, 1 - EPSILON)
    y_logit = np.log(y_clipped / (1 - y_clipped))

    # 2. Fit GAM
    gam = LinearGAM(s(0, n_splines=N_SPLINES, lam=LAMBDA)).fit(X_age, y_logit)

    # 3. Predict
    predictions_logit = gam.predict(age_range_pred)
    predictions_fraction = inverse_logit(predictions_logit) 

    # 4. Plot
    ax.plot(
        age_range_pred, 
        predictions_fraction,
        color=line_color, 
        linewidth=1.5, # Reduced line thickness to match smaller figure scale
        label=cell_type
    )

    # Inherits axes.titlesize (8pt) and font.sans-serif (Arial) automatically
    ax.set_title(f'{cell_type}', weight='bold') 
    ax.set_xlabel(f'Age ({AGE_COL})') # Inherits axes.labelsize (7pt)
    ax.set_ylabel('Fraction')           # Inherits axes.labelsize (7pt)
    ax.set_ylim(0, None)
    ax.legend(loc='upper right', frameon=False)

# Inherits global styling, kept slightly larger for main title
plt.suptitle('Cell Type Fractions Smoothed by Age using Linear GAM', weight='bold', y=1.03)
plt.tight_layout()

plt.savefig(
    '/home/levinsj/MultiOme/Analysis/pythonNotebooks/pythonPlots/cellFractionGAMs_overtime.pdf', 
    format='pdf', 
    bbox_inches='tight'
)
plt.show()
plt.close()

