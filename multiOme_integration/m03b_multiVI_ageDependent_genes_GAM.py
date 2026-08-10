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

from pygam import GAM, s, te, PoissonGAM, GammaGAM
from pygam.utils import OptimizationError


import seaborn as sns
import scipy


# In[2]:


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

sc.settings.figdir = "/home/levinsj/MultiOme/Analysis/pythonNotebooks/pythonPlots/GEX_GAM"


# ### scanpy variables

# In[3]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2
res = 2.0
n_features = 250000


# In[4]:


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

# In[5]:


sc.settings.verbosity = 3
sc.logging.print_header()
mpl.rcParams['figure.dpi'] = 450
min_distUMAP = 0.2


# In[6]:


inputSCVIfile = "/home/levinsj/MultiOme/MergedObjects/merged_postMultiVI_annotated_onlyGEX_human.h5ad"
output_dir = "/home/levinsj/MultiOme/Analysis/DEGs/"


# In[7]:


adata = sc.read_h5ad(inputSCVIfile)
adata.X = adata.layers["counts"]
adata.obs["age"] = adata.obs["sample_final"]
adata.obs["gestational_age"] = adata.obs["sample_final"]

adata.obs["age"] = adata.obs["age"].map(age_dict).astype('category')
values_to_replace = ['45/61', '49/68', '20/48']
adata.obs["age"] = adata.obs["age"].replace(values_to_replace, np.nan)

non_nan_age_mask = adata.obs['age'].notna()
adata_cleaned = adata[non_nan_age_mask, :].copy() 
adata_cleaned.obs['age'] = adata_cleaned.obs['age'].astype(float)

adata_cleaned.obs["gestational_age"] = adata_cleaned.obs["gestational_age"].map(gAge_dict).astype('category')
adata_cleaned.obs["gestational_age"] = adata_cleaned.obs["gestational_age"].str.slice(0, 2)
adata_cleaned.obs["gestational_age"] = adata_cleaned.obs["gestational_age"].replace("No", np.nan)
adata_cleaned.obs['gestational_age'] = adata_cleaned.obs['gestational_age'].astype(float)
adata_cleaned.obs['gestational_age'] = 0-((39-adata_cleaned.obs['gestational_age'])/52)

has_gestational_age = adata_cleaned.obs['gestational_age'].notna()

adata_cleaned.obs.loc[has_gestational_age, 'age'] = adata_cleaned.obs.loc[has_gestational_age, 'gestational_age']
age_data = adata_cleaned.obs['age'].dropna()


# #### Get cell types to test, so we can loop through them

# In[8]:


def plot_gene_expression_gam(goi: str, adata, age_max: float, cell_types: list = None, use_symlog_scale: bool = False, show_averages: bool = False, ignore_g_age: bool = False):
    """
    Plots GAM-fitted gene expression over 'age' for different cell types,
    including 95% confidence intervals and average expression per sample at each age.

    Args:
        goi (str): The gene of interest (must be in adata.var.index).
        adata: Your AnnData object containing expression and observation data.
               It should have 'cellType' and 'age' in adata.obs
               and gene names in adata.var.index.
        age_max (float): The maximum age to consider for plotting.
        cell_types (list, optional): A list of cell types to include in the plot.
                                     If None, defaults to ['PT', 'LOH', 'Endothelium', 'Podocyte', 'Stroma'].

    Raises:
        ValueError: If the gene of interest is not found in adata.var.index.
    """

    # --- GAM Model Parameters ---
    n_splines = 6 # Number of splines for the GAM model (can increase for more responsiveness)
    lam = 2     # Regularization parameter (decrease for more responsiveness)
    spline_order = 3 # Order of the spline interpolation

    if cell_types is None:
        cell_types = ['PT', 'LOH', 'DCT_CNT', 'Endothelium', 'Stroma']

    x_param = 'age' # The observation column to use as the independent variable (age)
    df = pd.DataFrame() # DataFrame to store GAM fit results for each cell type

    # --- Gene of Interest Validation ---
    if goi not in adata.var.index:
        raise ValueError(f"Gene '{goi}' not found in adata.var.index. Please check your gene name or adata.var.")
    goi_idx = adata.var.index.get_loc(goi) # Get the numerical index of the gene of interest

    # --- Plot Setup ---
    plt.figure(figsize=(6, 4)) # Create a new figure with a specified size for better visibility

    # Subset adata once for ages below age_max to avoid redundant operations

    if ignore_g_age:
        adata_filtered_age = adata[adata.obs["age"] > 0]
        print("ignoring gestational age")

    else:
        adata_filtered_age = adata
        print("using gestational age")

    # --- Iterate through Cell Types ---
    for cell_type in cell_types:
        # Create a subset for the current cell type
        adata_subset = adata_filtered_age[adata_filtered_age.obs['cellType'] == cell_type, :].copy()

        if adata_subset.n_obs == 0:
            print(f"Skipping {cell_type}: No observations found for this cell type.")
            continue

        # Extract 'age' and gene expression data for the current subset
        X_age = adata_subset.obs[[x_param]].to_numpy()
        y_expression = adata_subset.X[:, goi_idx].todense()
        y_expression[y_expression < 0] = 0 # Ensure expression values are non-negative for Poisson GAM
        y_expression = np.asarray(y_expression).flatten()

        # --- Calculate Average Expression per Unique Age ---
        # Create a temporary DataFrame to group by age and calculate mean expression
        temp_sample_df = pd.DataFrame({'age': X_age.flatten(), 'expression': y_expression})
        # Group by age and compute the mean expression for each unique age point
        avg_expression_per_age = temp_sample_df.groupby('age')['expression'].mean().reset_index()

        # --- Fit GAM Model for Current Cell Type ---
        current_model = PoissonGAM(s(0, spline_order=spline_order, n_splines=n_splines, lam = lam),
                                   max_iter=2000, fit_intercept=True, tol=0.001, verbose=False)
        current_model.fit(X_age, y_expression)

        # Generate points for plotting the GAM curve
        XX1 = np.linspace(-0.5, age_max+1, 500) # Create a smooth range of ages for prediction
        y_pred1 = current_model.predict(XX1) # Predict gene expression using the GAM model

        # Calculate 95% confidence intervals for the GAM prediction
        ci = current_model.confidence_intervals(XX1, width=0.95)

        # Store results in the DataFrame for potential external use (optional)
        df[f'GEX_{cell_type}'] = pd.Series(y_pred1.flatten(), name=f'GEX_{cell_type}')
        df[f'age_{cell_type}'] = pd.Series(XX1, name=f'age_{cell_type}')
        df[f'GEX_{cell_type}_lower'] = pd.Series(ci[:, 0].flatten(), name=f'GEX_{cell_type}_lower')
        df[f'GEX_{cell_type}_upper'] = pd.Series(ci[:, 1].flatten(), name=f'GEX_{cell_type}_upper')

        line, = plt.plot(df[f'age_{cell_type}'], df[f'GEX_{cell_type}'], label=cell_type)
        line_color = line.get_color()

        plt.fill_between(df[f'age_{cell_type}'], df[f'GEX_{cell_type}_lower'], df[f'GEX_{cell_type}_upper'], alpha=0.2, color=line_color)
        if show_averages:
            plt.scatter(avg_expression_per_age['age'], avg_expression_per_age['expression'], color=line_color, s=20, alpha=0.6, marker='o', edgecolor='none', label=f'{cell_type} Avg. Sample GEX')


    if df.empty:
        print("No data was processed for plotting. Please check your inputs and cell types.")
        plt.close() # Close any empty plot that might have been created
        return df # Return the empty DataFrame

    plt.title(f'GAM Fit for {goi} by Cell Type with 95% Confidence Intervals and Sample Averages')
    plt.xlabel('Age')
    plt.ylabel(f'{goi} Expression')
    if use_symlog_scale:
        plt.xscale('symlog', linthresh=5)
    if ignore_g_age:
        plt.xlim(0, age_max)
    else:
        plt.xlim(-0.6, age_max)
        plt.axvline(x=0, color='black', linestyle='--', linewidth=1.0, label='Birth') # Vertical line at age 0 for birth

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.grid(True, linestyle='-', alpha=0.7)
    plt.suptitle(f'{goi} ChromVAR GAM', fontsize=16, weight='bold', y=1.02)

    fig = plt.gcf()
    fig.tight_layout()
    fig.savefig(f'/home/levinsj/MultiOme/Analysis/pythonNotebooks/pythonPlots/GEX_GAM/GEX_GAM_{goi}_GAM_ages_0_{age_max}.pdf', format='pdf', bbox_inches='tight')
    plt.close(fig)
    del df
    del temp_sample_df


# ### Fetal Specific genes

# In[9]:


goi = ["APOL1","PLA2R1","NPHS1","THSD7A","NELL1","NPHS2"]

for i in goi: 

    try:
        plot_gene_expression_gam(i, adata_cleaned, 80, use_symlog_scale=False, show_averages=False, ignore_g_age = False, cell_types = ["Podocyte","Endothelium","Stroma"])

    except Exception as e:
        print(f"An error occurred during the GAM analysis for {i}")
        print(e)


# In[10]:


goi = ["NPHS1", "NPHS2", "PODXL"]

for i in goi: 

    try:
        plot_gene_expression_gam(i, adata_cleaned, 20, use_symlog_scale=False, show_averages=False, ignore_g_age = False, cell_types = ["Podocyte","Endothelium","Stroma"])

    except Exception as e:
        print(f"An error occurred during the GAM analysis for {i}")
        print(e)


# In[11]:


### Fetal Specific Genes


# In[12]:


goi = ["HIF3A","HMGA2"]

for i in goi: 

    try:
        plot_gene_expression_gam(i, adata_cleaned, 3, use_symlog_scale=False, show_averages=False, ignore_g_age = False, cell_types = ["PT"])

    except Exception as e:
        print(f"An error occurred during the GAM analysis for {i}")
        print(e)


# In[13]:


goi = ["GSTA2","HNF4A"]

for i in goi: 

    try:
        plot_gene_expression_gam(i, adata_cleaned, 80, use_symlog_scale=False, show_averages=False, ignore_g_age = False, cell_types = ["PT"])

    except Exception as e:
        print(f"An error occurred during the GAM analysis for {i}")
        print(e)


# goi = ["HMGA2","WT1","SOX4","SOX11", "SOX12", "EYA1", "HIF3A", "VCAN"]
# 
# for i in goi: 
# 
#     try:
#         plot_gene_expression_gam(i, adata_cleaned, 3, use_symlog_scale=False, show_averages=False, ignore_g_age = False)
# 
#     except Exception as e:
#         print(f"An error occurred during the GAM analysis for {i}")
#         print(e)

# goi = ["ESR1", "KL","HNF4A","HNF1A","TMEM63A","SLC17A9", "HNF4G", "HNF1B", "GSTA1","GSTA2","HIF3A","VCAN","PCDH15","RDX","VCAM1"]
# for i in goi: 
#     try:
#         plot_gene_expression_gam(i, adata_cleaned, 80, use_symlog_scale=False, show_averages=False, ignore_g_age = False)
# 
#     except Exception as e:
#         print(f"An error occurred during the GAM analysis for {i}")
#         print(e)

# goi = ["IRX3","ANKFN", "CDYL2", "MBOAT1", "KCND3", "NCNK2", "ABDG2"]
# for i in goi: 
#     try:
#         plot_gene_expression_gam(i, adata_cleaned, 80, use_symlog_scale=False, show_averages=False, ignore_g_age = False)
#         plot_gene_expression_gam(i, adata_cleaned, 3, use_symlog_scale=False, show_averages=False, ignore_g_age = False)
# 
#     except Exception as e:
#         print(f"An error occurred during the GAM analysis for {i}")
#         print(e)

# ### Adult Tx changes

# goi = ["MEIS2","KLF9","ARNT2","EYA2","FOXO1","PLCG2","KL","LDB2"]
# 
# for i in goi: 
#     try:
#         plot_gene_expression_gam(i, adata_cleaned, 80, use_symlog_scale=False, show_averages=False, ignore_g_age = False)
# 
#     except Exception as e:
#         print(f"An error occurred during the GAM analysis for {i}")
#         print(e)
