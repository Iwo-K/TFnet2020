# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.1.6
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # TFnet - plotting gene expression of gene modules

# %%
import numpy as np
import pandas as pd
import scanpy.api as sc
import matplotlib.pyplot as plt
import tables as tb
import scipy as scipy
import dotscore

#Specifying random seed
import random
random.seed(123)

sc.set_figure_params(color_map = 'viridis', dpi_save = 350)
sc.settings.verbosity = 3

# %%
#To ensure 1:1 aspect ratio
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (4, 4)

# %%
import pyfuns

# %%
sc.settings.figdir = './module_plot/'

# %% [markdown]
# ## Loading data

# %% [markdown]
# ### Loading the expression changes - DE data

# %%
dedf = pd.read_csv("./NET/NETdf_filt.csv", index_col=0) 

# %% [markdown]
# Extracting all observed changes in gene expression, which will use as a set for simulations to estimate DoT z-score

# %%
allfolds = dedf['log2FoldChange'].values

# %% [markdown]
# Loading genes expressed in Hoxb8-FL cells (also required for the simulation)

# %%
exprdata = pd.read_csv("./procdata/TFnet_counts_QC_expressedset.csv", index_col = 0)
allgenes = exprdata.index.values

# %% [markdown]
# Gene annotation

# %%
genedata = pd.read_csv("./data/TFnet_data/genedata.csv", index_col = 0)

# %% [markdown]
# ## Loading reference scRNA-Seq landscapes

# %% [markdown]
# ### LKdata

# %% [markdown]
# Loading pre-processed data from Dahlin et al. Blood 2018 paper. It contains cells from Lin-, Kit+ and Lin-, Kit+, Sca1+ populations. Thus the landscape encompasses stem cells as well as multi and uni-potent progenitors. We also set the desired UMAP coordinates.

# %%
LKadata = sc.read('./data/TFnet_data/LKLSKdata_QC_processed_small.h5ad')
LKadata.var.index = LKadata.var.geneid
LKadata.obsm['X_umap'] = LKadata.obsm['X_umap_noG']

# %%
sc.pp.scale(LKadata)

# %% [markdown]
# ### HCA BM

# %% [markdown]
# We load the pre-processed human bone marrow mononuclear cell data. Data was obtained from the Human Cell Atlas (https://data.humancellatlas.org/explore/projects/cc95ff89-2e68-4a08-a234-480eca21ce79). This dataset contains ~250,000 cells spanning entire trajectories from stem cell up to fully differentiated cells.
#
# **WARNING** this dataset is large and its processing requires >20GB RAM memory

# %%
HCAdata = sc.read('./data/TFnet_data/hcaBM_processed_small.h5ad')
HCAdata = HCAdata[:, ~(HCAdata.var.mouseid == 'None')].copy()
HCAdata.var.index = HCAdata.var.mouseid

# %%
sc.pp.scale(HCAdata)


# %% [markdown]
# ## Plotting expression of specific groups of genes

# %% [markdown]
# Convenience functions to extract co-regulated genes from the network data (dedf data frame)

# %%
def get_genes2(df, 
               which, 
               qs = [1], 
               direction = 'none', 
               mode = 'none', 
               plot = True, 
               save_genes = None, 
               genedata = genedata):
    
    coregs = df.loc[df.ko.isin(which),:]
    coregs = coregs.pivot(index = 'target', columns = 'ko', values = 'log2FoldChange')
    coregs = coregs.dropna()
    print("There are " + str(len(coregs.index)) + " genes in the overlap")

    if(plot):
        ax = plt.scatter(coregs[which[0]], coregs[which[1]], alpha = 0.7)
        plt.xlabel(which[0])
        plt.ylabel(which[1])
        plt.axvline(x=0)
        plt.axhline(y=0)
        plt.show()
    
    coregs['quadrants'] = 0
    for i in coregs.index:
        if coregs.loc[i,which[0]] > 0 and coregs.loc[i,which[1]] > 0:
            coregs.loc[i, 'quadrants'] = 1
        elif coregs.loc[i,which[0]] < 0 and coregs.loc[i,which[1]] < 0:
            coregs.loc[i, 'quadrants'] = 3
        elif coregs.loc[i,which[0]] > 0 and coregs.loc[i,which[1]] < 0:
            coregs.loc[i, 'quadrants'] = 4
        elif coregs.loc[i,which[0]] < 0 and coregs.loc[i,which[1]] > 0:
            coregs.loc[i, 'quadrants'] = 2
        else:
            print('There is something wrong with the values')
    

    
    out = coregs.index[coregs.quadrants.isin(qs)]
    print('Selected: ', str(len(out)), ' genes')
    
    if save_genes != None:
        coregs['symbol'] = [genedata.symbol[x] if x in genedata.index else None for x in coregs.index ]
        coregs.to_csv(save_genes)
    return(out)

# %% [markdown]
# ## Plotting expression target shared between pairs of TFs

# %% [markdown]
# We select genes co-regulated by Gata3 and Ebf1 and split them into 4 groups (quadrants) dendeing on the pattern of Gata3 and Ebf1 coregulation (+/+, +/-, -/+, -/-).
# We can plot the sum of scaled expression values to observe where these genes are 

# %%
pairs = [['Gata3_Plate5', 'Ebf1_Plate16']]
for i in pairs:
    fig, axes = plt.subplots(nrows=2, ncols=4)
    fig.set_size_inches(20, 10)
    for j in [1,2,3,4]:
        
        x = get_genes2(dedf, which = i, qs = [j], plot = False, save_genes = './module_plot/' + '__'.join(i) + '.csv')
        name = '_'.join(i) + '_q' + str(j) + '_' + str(len(x)) + 'genes'
        
        LKadata.obs[name] = pyfuns.get_module_Zscore(LKadata, x, allgenes = allgenes, simno = 500)
        LKadata.obs[name] = dotscore.qfilt(LKadata.obs[name], 0.999, 0.001)
        cmap0 = dotscore.cmap_RdBu(LKadata.obs[name])
        sc.pl.umap(LKadata, color = name, 
                   alpha = 0.7,
                   ax = axes[0,j-1], show = False, cmap = cmap0, frameon = False)
        
        HCAdata.obs[name] = pyfuns.get_module_Zscore(HCAdata, x, allgenes = allgenes, simno = 500)
        HCAdata.obs[name] = dotscore.qfilt(HCAdata.obs[name], 0.999, 0.001)
        cmap0 = dotscore.cmap_RdBu(HCAdata.obs[name])
        sc.pl.umap(HCAdata, color = name, alpha = 0.7, ax = axes[1,j-1], 
                   show = False, 
                   cmap = cmap0,
                   title = '', frameon = False)
    plt.savefig('./module_plot/' + '__'.join(i) + '.pdf')
