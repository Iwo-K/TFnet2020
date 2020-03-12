# -*- coding: utf-8 -*-
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
# # TFnet - calculating DoT scores for TF perturbations

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
# ## Loading reference scRNA-Seq landscapes

# %% [markdown]
# ### LKdata

# %% [markdown]
# Loading pre-processed data from Dahlin et al. Blood 2018 paper. It contains cells from Lin-, Kit+ and Lin-, Kit+, Sca1+ populations. Thus the landscape encompasses stem cells as well as multi and uni-potent progenitors. We also set the desired UMAP coordinates.

# %%
LKadata = sc.read('./data/TFnet_data/LKLSKdata_QC_processed_small.h5ad')
LKadata.var.index = LKadata.var.geneid
LKadata.obsm['X_umap'] = LKadata.obsm['X_umap_noG']

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

# %% [markdown]
# ### DoT score - LKdata

# %% [markdown]
# Scaling to the cluster most closely resembling Hoxb8-FL cells - cluster 4 (see previous script for details)

# %%
LKadata2 = LKadata.copy()
means4 = LKadata2[LKadata2.obs.leiden_noG == '4',:].X.mean(axis = 0)
dotscore.custom_scale(LKadata2, mean = means4)

# %% [markdown]
# Computing DoT scores for each perturbation

# %%
#Setting the output folder
sc.settings.figdir = './TFnet_DoTscore/LKprojections_CL4scaled/'
kos = list(set(dedf['ko'].tolist()))
for i in kos:
    deI = dedf[dedf['ko'] == i]
    scorename = i + 'score'
    LKadata2.obs[scorename + "Z"] = dotscore.get_DoTscore(LKadata2, 
                                                 de = deI, 
                                                 allgenes = allgenes, 
                                                 allfolds = allfolds,
                                                 zscore = True)
    LKadata2.obs[scorename + "Z_qfilt"] = dotscore.qfilt(LKadata2.obs[scorename + "Z"], 0.999, 0.001)
    
    c1 = dotscore.cmap_RdBu(LKadata2.obs[scorename + "Z_qfilt"])

    filename = 'LKdata_' + scorename + "DoTscore" + '.pdf'
    sc.pl.umap(LKadata2, color=[scorename + "Z_qfilt"], show = False, alpha = 0.7, color_map = c1, frameon = False, save = filename)
    
    leiden_scores = dotscore.get_genescore_pergroup(LKadata2, deI, group = 'leiden_noG', sortby = '0', gene_symbols = 'symbol')
    leiden_scores = leiden_scores.loc[leiden_scores.sum(axis = 1) != 0, :]
    leiden_scores.to_csv('./TFnet_DoTscore/LKprojections_CL4scaled/' + 'LKdata_' + scorename + "_DoTscore_per_leidennoG" + '.csv')

# %% [markdown]
# ### HCA FC scores for each using (scaled to the CD34+ compartment)

# %% [markdown]
# In case of human cells, we will use the cluster number 22, which contains the CD34+ (here using mouse id: ENSMUSG00000016494) cells as a point of origin.

# %%
sc.pl.umap(HCAdata, color = ['louvain', 'ENSMUSG00000016494'], legend_loc = 'on data')

# %%
means22 = HCAdata[HCAdata.obs.louvain == '22',:].X.mean(axis = 0)
dotscore.custom_scale(HCAdata, mean = means22)

# %% [markdown]
# Computing DoT scores for each perturbation

# %%
#Setting the output folder
sc.settings.figdir = './TFnet_DoTscore/HCAprojections_CL22scaled/'
kos = list(set(dedf['ko'].tolist()))
for i in kos:
    deI = dedf[dedf['ko'] == i]
    scorename = i + 'score'
    HCAdata.obs[scorename + "Z"] = dotscore.get_DoTscore(HCAdata, 
                                                 de = deI, 
                                                 allgenes = allgenes, 
                                                 allfolds = allfolds,
                                                 zscore = True)
    HCAdata.obs[scorename + "Z_qfilt"] = dotscore.qfilt(HCAdata.obs[scorename + "Z"], 0.999, 0.001)
    
    c1 = dotscore.cmap_RdBu(HCAdata.obs[scorename + "Z_qfilt"])
    
    filename = 'HCdata_' + scorename + "DoTscore" + '.pdf'
    sc.pl.umap(HCAdata, color=[scorename + "Z_qfilt"], show = False, alpha = 0.7, color_map = c1, frameon = False, save = filename)

    louvain_scores = dotscore.get_genescore_pergroup(HCAdata, deI, group = 'louvain', sortby = '0', gene_symbols = 'mousesymbol')
    louvain_scores = louvain_scores.loc[louvain_scores.sum(axis = 1) != 0, :]
    louvain_scores.to_csv('./TFnet_DoTscore/HCAprojections_CL22scaled/' + 'HCAdata_' + scorename + "_DoTscore_per_louvain" + '.csv')
