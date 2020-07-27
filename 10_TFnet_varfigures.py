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
# # TFnet - landscape figures

# %%
import numpy as np
import pandas as pd
import scanpy.api as sc
import matplotlib.pyplot as plt
import tables as tb
import scipy as scipy
import dotscore

sc.set_figure_params(color_map = 'viridis', dpi_save = 350)
sc.settings.verbosity = 3

# %%
#To ensure 1:1 aspect ratio
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (4, 4)

# %%
sc.settings.figdir = './varfigures/'


# %%
def annotate_clusters(adataobj, df, cluster_column, anno_column, color_column):
    adataobj.obs[anno_column] = [df.loc[df.index[int(x)], anno_column] for x in adataobj.obs[cluster_column]]
    adataobj.obs['anno'] = adataobj.obs['anno'].astype('category')
    
    d = {df.loc[df.index[x], anno_column] : df.loc[df.index[x],color_column] for x in range(len(df.index))}
    cats = adataobj.obs[anno_column].cat.categories.to_list()
    adataobj.uns[anno_column + '_colors'] = [d[x] for x in cats]
    


# %% [markdown]
# ## LK data

# %% [markdown]
# Loading pre-processed data from Dahlin et al. Blood 2018 paper. It contains cells from Lin-, Kit+ and Lin-, Kit+, Sca1+ populations. Thus the landscape encompasses stem cells as well as multi and uni-potent progenitors

# %%
LKadata = sc.read('./data/TFnet_data/LKLSKdata_QC_processed_small.h5ad')
LKadata.obsm['X_umap'] = LKadata.obsm['X_umap_noG']

# %% [markdown]
# Scaling to the cluster 1

# %%
means1 = LKadata[LKadata.obs.leiden_noG == '1',:].X.mean(axis = 0)
dotscore.custom_scale(LKadata, mean = means1)

# %% [markdown]
# ### Plotting DoT score example

# %%
LKadata.obs['Klf1_expr'] = LKadata.X[:, LKadata.var.index == 'Klf1']
cmap0 = dotscore.cmap_RdBu(LKadata.obs['Klf1_expr'])
sc.pl.umap(LKadata, color = 'Klf1_expr', cmap = cmap0, use_raw = False, save = '_LKdata_Klf1expr.pdf', title = '', frameon = False)

LKadata.obs['Klf1_expr_pos1_7'] = LKadata.X[:, LKadata.var.index == 'Klf1']*1.7
cmap0 = dotscore.cmap_RdBu(LKadata.obs['Klf1_expr_pos1_7'])
sc.pl.umap(LKadata, color = 'Klf1_expr_pos1_7', cmap = cmap0, use_raw = False, save = '_LKdata_Klf1expr_pos1_7.pdf', title = '', frameon = False)


LKadata.obs['Elane_expr'] = LKadata.X[:, LKadata.var.index == 'Elane']
cmap0 = dotscore.cmap_RdBu(LKadata.obs['Elane_expr'])
sc.pl.umap(LKadata, color = 'Elane_expr', cmap = cmap0, use_raw = False, save = '_LKdata_Elaneexpr.pdf', title = '', frameon = False)

LKadata.obs['Elane_expr_neg0_5'] = LKadata.X[:, LKadata.var.index == 'Elane']*-0.5
cmap0 = dotscore.cmap_RdBu(LKadata.obs['Elane_expr_neg0_5'])
sc.pl.umap(LKadata, color = 'Elane_expr_neg0_5', cmap = cmap0, use_raw = False, save = '_LKdata_Elaneexpr_neg0_5.pdf', title = '', frameon = False)


LKadata.obs['Dntt_expr'] = LKadata.X[:, LKadata.var.index == 'Dntt']
cmap0 = dotscore.cmap_RdBu(LKadata.obs['Dntt_expr'])
sc.pl.umap(LKadata, color = 'Dntt_expr', cmap = cmap0, use_raw = False, save = '_LKdata_Dnttexpr.pdf', title = '', frameon = False)

LKadata.obs['Dntt_expr_neg1_25'] = LKadata.X[:, LKadata.var.index == 'Dntt']*-1.25
cmap0 = dotscore.cmap_RdBu(LKadata.obs['Dntt_expr_neg1_25'])
sc.pl.umap(LKadata, color = 'Dntt_expr_neg1_25', cmap = cmap0, use_raw = False, save = '_LKdata_Dnttexpr_neg1_25.pdf', title = '', frameon = False)

LKadata.obs['sum_expr'] = LKadata.obs['Dntt_expr']*-1.25 + LKadata.obs['Elane_expr']*-0.5 + LKadata.obs['Klf1_expr']*1
cmap0 = dotscore.cmap_RdBu(LKadata.obs['sum_expr'])
sc.pl.umap(LKadata, color = 'sum_expr', cmap = cmap0, use_raw = False, save = '_LKdata_sumexpr.pdf', title = '', frameon = False)

# %% [markdown]
# ## Plotting landscape annotation

# %% [markdown]
# ### LKdata

# %% [markdown]
# Loading annotation

# %%
sc.pl.umap(LKadata, color = 'leiden_noG', legend_loc = 'on data')

# %%
xrecol = pd.read_csv('./data/TFnet_data/LKdata_clustercolors_recolored.csv')

# %%
annotate_clusters(LKadata, xrecol, 'leiden_noG', anno_column = 'anno', color_column = 'recolored6')
sc.pl.umap(LKadata, color = 'anno', legend_loc = None, frameon = False, title = '', save = '_LKLSK_clusters_recolored6', alpha = 0.7, size = 8)

# %% [markdown]
# ## HCA data

# %% [markdown]
# Loading data and annotation

# %%
HCAdata = sc.read('./data/TFnet_data/hcaBM_processed_small.h5ad')

# %%
xrecol = pd.read_csv('./data/TFnet_data/HCAdata_clustercolors_recolored.csv')

# %%
annotate_clusters(HCAdata, xrecol, 'leiden', anno_column = 'anno', color_column = 'recolored4')
sc.pl.umap(HCAdata, color = 'anno', legend_loc = None, frameon = False, title = '', save = '_HCA_clusters_recolored4', alpha = 0.7, size = 2)

# %%
