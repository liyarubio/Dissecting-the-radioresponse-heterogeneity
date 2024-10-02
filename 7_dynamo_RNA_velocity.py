# %%
import dynamo as dyn
import os
import sys
import scanpy as sc
import scvelo as scv

# %%
dyn.get_all_dependencies_version()

# %%
dyn.configuration.set_figure_params('dynamo', background='white')

# %%
#adata2 = dyn.sample_data.zebrafish()

# %%
#adata = sc.read("./Ctrl.bm2.TSNE.h5ad")
adata = sc.read("/media/liyaru/LYR1/TH/velocyto/NC.h5ad")

# %%
adata

# %%
dyn.pp.recipe_monocle(adata)

# %%
dyn.tl.dynamics(adata, model='stochastic', cores=3)

# %%
dyn.tl.reduceDimension(adata)
dyn.pl.tsne(adata, color='celltype')

# %%
dyn.pl.umap(adata, color='celltype')

# %%
# dyn.tl.gene_wise_confidence(adata, group='celltype', lineage_dict={'HSPC-1': ['NP']})

# %%
dyn.pl.phase_portraits(adata, genes=adata.var_names[adata.var.use_for_dynamics][:4], figsize=(6, 4), color='celltype')

# %%
dyn.tl.cell_velocities(adata, method='pearson', other_kernels_dict={'transform': 'sqrt'})

# %%
dyn.tl.cell_wise_confidence(adata)
# dyn.tl.confident_cell_velocities(adata, group='group', lineage_dict={'Progenitor': ['terminal_cell_state']},)
# dyn.tl.confident_cell_velocities(adata, group='celltype', lineage_dict={'HSPC-1': ['NP','cDC',]})

# %%
Color_key= {'G1':'#FFFFB3',
           'S':'#8DD3C7',
           'G2M':'#BEBADA'}

# %%
dyn.pl.streamline_plot(adata, color=['celltype'], basis='umap', show_legend='on data', show_arrowed_spines=True,
                       color_key = Color_key,
                       #show_legend='right',
                      figsize = [8,8])

# %%
dyn.pl.streamline_plot(adata, color=['celltype'], basis='umap',  show_arrowed_spines=True,
                      color_key = Color_key,
                      save_show_or_return="save",
                      save_kwargs = {"path":"NC_UMAP",
                                     "ext":"png",
                                    "dpi":500},
                      show_legend = None,
                       label=None,
                      figsize = [6,6])

# %%
# test save
# plt.semilogy([3,1,4,-1,5,9])
# plt.savefig('semilogy.pdf')

# %%
dyn.vf.VectorField(adata, basis='umap', M=1000, pot_curl_div=True)

# %%
dyn.pl.topography(adata, basis='umap', background='white', 
                  color=['ntr','celltype'], streamline_color='black', 
                  show_legend='on data', frontier=True)

# %%
dyn.pl.umap(adata,  color='umap_ddhodge_potential', frontier=True)

# %%
dyn.pl.umap(adata,  color='umap_ddhodge_potential', frontier=True,
           save_show_or_return="save",
           save_kwargs = {"path":"UMAP",
                          "ext":"pdf",
                          "dpi":200},
                      figsize = [8,8]) 

# %%
stem_marker = ["CCNB1"]

# %%
import numpy as np
dyn.pl.scatters(adata, x=np.repeat('umap_ddhodge_potential', 9), pointsize=0.25, alpha=0.8, 
                y=stem_marker, 
                layer='X_spliced', color='celltype',
                ncols=3, background='white', 
                #color_key = Color_key,
                figsize=(9, 4))

# %%
dyn.tl.cell_velocities(adata, basis='pca')
dyn.vf.VectorField(adata, basis='pca')
dyn.vf.speed(adata, basis='pca')
dyn.vf.curl(adata, basis='umap')
dyn.vf.divergence(adata, basis='pca')
dyn.vf.acceleration(adata, basis='pca')
dyn.vf.curvature(adata, basis='pca')

# %%
import matplotlib.pyplot as plt

fig1, f1_axes = plt.subplots(ncols=2, nrows=2, constrained_layout=True, figsize=(12, 8))
f1_axes
f1_axes[0, 0] = dyn.pl.cell_wise_vectors(adata, color='speed_pca', pointsize=0.5, alpha = 0.7, ax=f1_axes[0, 0], quiver_length=6, quiver_size=6, save_show_or_return='return')
f1_axes[0, 1] = dyn.pl.grid_vectors(adata, color='divergence_pca', ax=f1_axes[0, 1], quiver_length=12, quiver_size=12, save_show_or_return='return')
f1_axes[1, 0] = dyn.pl.streamline_plot(adata, color='acceleration_pca', ax=f1_axes[1, 0], save_show_or_return='return')
f1_axes[1, 1] = dyn.pl.streamline_plot(adata, color='curvature_pca', ax=f1_axes[1, 1], save_show_or_return='return')
plt.show()

# %%
dyn.configuration.set_figure_params('dynamo', background='black')

# %%
fig1, f1_axes = plt.subplots(ncols=2, nrows=2, constrained_layout=True, figsize=(12, 8))
f1_axes
f1_axes[0, 0] = dyn.pl.cell_wise_vectors(adata, color='speed_pca', pointsize=0.1, alpha = 0.7, ax=f1_axes[0, 0], quiver_length=6, quiver_size=6, save_show_or_return='return', background='black')
f1_axes[0, 1] = dyn.pl.grid_vectors(adata, color='divergence_pca', ax=f1_axes[0, 1], quiver_length=12, quiver_size=12, save_show_or_return='return', background='black')
f1_axes[1, 0] = dyn.pl.streamline_plot(adata, color='acceleration_pca', ax=f1_axes[1, 0], save_show_or_return='return', background='black')
f1_axes[1, 1] = dyn.pl.streamline_plot(adata, color='curvature_pca', ax=f1_axes[1, 1], save_show_or_return='return', background='black')
plt.show()

# %%
#progenitor = adata.obs_names[adata.obs.celltype.isin(['G2M'])]
progenitor = adata.obs_names[adata.obs.celltype.isin(['G1'])]
len(progenitor)

# %%
dyn.pd.fate(adata, basis='umap', init_cells=progenitor, interpolation_num=100,  direction='forward',
   inverse_transform=False, average=False, cores=3)

# %%
%%capture
fig, ax = plt.subplots()
ax = dyn.pl.topography(adata, color='celltype', ax=ax, save_show_or_return='return')

# %%
%%capture
instance = dyn.mv.StreamFuncAnim(adata=adata, color='celltype', ax=ax)

# %%
import matplotlib
matplotlib.rcParams['animation.embed_limit'] = 2**128 # Ensure all frames will be embedded.

from matplotlib import animation
import numpy as np

anim = animation.FuncAnimation(instance.fig, instance.update, init_func=instance.init_background,
                               frames=np.arange(100), interval=100, blit=True)
from IPython.core.display import display, HTML
HTML(anim.to_jshtml()) # embedding to jupyter notebook.

# %%
dyn.mv.animate_fates(adata, color='celltype', basis='umap', n_steps=100, fig=fig, ax=ax,
                     save_show_or_return='save', logspace=True, max_time=None)

# %%
%%capture
fig, ax = plt.subplots()
ax = dyn.pl.topography(adata, color='celltype', ax=ax, save_show_or_return='return')
dyn.mv.animate_fates(adata, color='celltype', basis='umap', n_steps=100, fig=fig, ax=ax,
                     save_show_or_return='save', logspace=True, max_time=None)


