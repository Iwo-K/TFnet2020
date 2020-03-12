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
# # TFnet - projecting Hoxb8-FL cells onto scRNAseq landscapes

# %%
import numpy as np
import pandas as pd
import scanpy.api as sc
import matplotlib.pyplot as plt
import tables as tb
import scipy as scipy

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
# Importing the 

# %%
import pyfuns as pyfuns

# %% [markdown]
# Output files directory

# %%
sc.settings.figdir = './CELLproj/'

# %% [markdown]
# ## Loading data

# %% [markdown]
# ### scRNA-Seq landscape

# %% [markdown]
# Loading pre-processed data from Dahlin et al. Blood 2018 paper. It contains cells from Lin-, Kit+ and Lin-, Kit+, Sca1+ populations. Thus the landscape encompasses stem cells as well as multi and uni-potent progenitors

# %%
LKdata_raw = sc.read('./data/TFnet_data/LKLSKdata_QC_raw.h5ad')
LKdata = sc.read('./data/TFnet_data/LKLSKdata_QC_processed_small.h5ad')

# %%
LKdata_raw

# %%
LKdata

# %% [markdown]
# Replacing the nlog scaled values with raw counts - these will be normalised and scaled together with the Hoxb8-FL data.

# %%
LKdata_raw = LKdata_raw[:, LKdata.var.index].copy()
LKdata.X = LKdata_raw.X
LKdata.var.index = LKdata.var.geneid

# %% [markdown]
# Removing cell cycle genes (this step is also performed to obtain the visualisation of scRNA-Seq data)

# %%
LKdata = LKdata[:, ~LKdata.var.isCC]

# %% [markdown]
# Loading the desired UMAP coordinates

# %%
LKdata.obsm['X_umap'] = LKdata.obsm['X_umap_noG']

# %% [markdown]
# ### Loding Hoxb8-FL scRNA-Seq data

# %% [markdown]
# We are loading scRNA-Seq data obtained for the WT Hoxb8-FL cells from Basilico et al. Nat Comm. 2020 paper.

# %%
Hoxb8 = sc.read("data/TFnet_data/Hoxb8_FL_WT_counts_QC.csv", cache = False).transpose()

# %% [markdown]
# ## Pre-processing data

# %% [markdown]
# Data is log-normalised and pre-transformed together.

# %%
target, ref = pyfuns.process_together(Hoxb8, LKdata, use_vargenes = True)

# %% [markdown]
# We find the nearest neighbors between the dataset in the PCA space calculated on the reference datasets. Subsequently we label the cells on the reference landscape (LK data) according to the number of neighbours found in the target data (here the Hoxb8-FL data).

# %%
LKdata.obs['Hoxb8_proj'] = pyfuns.project_to_ref(target, ref, k = 5)

# %% [markdown]
# Plotting the neighbours

# %%
from matplotlib.colors import LinearSegmentedColormap
cmap2 = LinearSegmentedColormap.from_list('mycmap', [(0, '#DDD3D3'),
                                                        (0.001, 'blue'),
                                                        (1, 'red')])

# %%
sc.pl.umap(LKdata, color = ['Hoxb8_proj'], cmap = cmap2, save = '_Hoxb8_proj_ontoLK_betterscale.pdf', frameon = False, size = 24)

# %%
from matplotlib.colors import LinearSegmentedColormap
cmap2 = LinearSegmentedColormap.from_list('mycmap', [(0, '#DDD3D3'),
                                                        (0.001, 'tan'),
                                                        (1, 'blue')])

# %%
sc.pl.umap(LKdata, color = ['Hoxb8_proj'], cmap = cmap2, save = '_Hoxb8_proj_ontoLK_betterscale.pdf', frameon = False, size = 64)

# %% [markdown]
# Plotting the clusters identified in the reference data

# %%
sc.pl.umap(LKdata, color = 'leiden_noG', legend_loc = 'on data')

# %% [markdown]
# Neighbouring cells span cluster 15 and cluster 4. As cluster 15 is the terminal cluster, cluster 4 seems to be suitable population reflecting the Hoxb8-FL cells.
