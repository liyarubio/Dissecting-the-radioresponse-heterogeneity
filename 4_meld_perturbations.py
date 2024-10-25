
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import graphtools as gt
import phate
import magic
import scprep
import meld
import cmocean
import sklearn
import scipy
import seaborn as sns

# setting defaults for matplotlib font sizes
import matplotlib.pyplot as plt
plt.rc('font', size=14)

# making sure plots & clusters are reproducible
np.random.seed(42)

get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import diffxpy.api as de


# In[2]:


data = pd.read_csv("/home/liyaru/res/expr_count.csv")


# In[3]:


data
#row gene column cell


# In[4]:


data = data.T


# In[5]:


data = data.T


# In[6]:


data= data.drop(['Unnamed: 0'])


# In[7]:


data
#row gene column cell


# In[8]:


metadata = pd.read_csv("/home/liyaru/res/meta_data.csv",index_col=0)


# In[9]:


metadata
#row cell


# In[10]:


#data = scprep.filter.filter_rare_genes(data)


# In[11]:


#scprep.plot.plot_library_size(data, cutoff=15000)


# In[12]:


#data, metadata = scprep.filter.filter_library_size(
#    data, metadata, cutoff=15000, 
#    keep_cells='below')


# In[13]:


#scprep.plot.plot_gene_set_expression(data, genes=['LOC101885394'], log='y', cutoff=164)


# In[14]:


data_libnorm, libsize = scprep.normalize.library_size_normalize(data, return_library_size=True)


# In[15]:


metadata['library_size'] = libsize
metadata.head()


# In[16]:


data_sqrt = np.sqrt(data_libnorm)


# In[17]:


fig, ax = plt.subplots(1)

groups, counts = np.unique(metadata['orig.ident'], return_counts=True)
for i, c in enumerate(counts):
    ax.bar(i, c)
    
ax.set_xticks(np.arange(i+1))
ax.set_xticklabels(groups)
ax.set_ylabel('# cells')

fig.tight_layout()


# In[18]:


data_pca = scprep.reduce.pca(data_sqrt)

phate_op = phate.PHATE(knn=10, decay=10, n_jobs=-1)
data_phate = phate_op.fit_transform(data_pca)


# In[19]:


scprep.plot.scatter2d(data_phate, c=metadata['orig.ident'], 
                      legend_anchor=(1,1), figsize=(6,5), s=10, label_prefix='PHATE', ticks=False)


# In[20]:


scprep.plot.scatter2d(data_phate, c=metadata['comb.cluster'], 
                      legend_anchor=(1,1), figsize=(6,5), s=10, label_prefix='PHATE', ticks=False)


# In[21]:


scprep.plot.scatter2d(data_phate, c=metadata['Phase'], 
                      legend_anchor=(1,1), figsize=(6,5), s=10, label_prefix='PHATE', ticks=False)


# In[22]:


metadata.columns


# In[23]:


metadata['sample_phase']


# In[24]:


metadata['Phase']


# In[25]:


scprep.plot.scatter2d(data_phate, c=metadata['orig.ident'], cmap=cmocean.cm.phase, 
                      legend_anchor=(1,1), figsize=(5,5), s=10, label_prefix='PHATE', ticks=False)


# In[26]:


metadata['IR'] = [1 if sl.startswith('IR') else 0 for sl in metadata['orig.ident']]


# In[27]:


metadata['sample']=metadata['orig.ident'] 


# In[28]:


metadata['sample2']=metadata['sample']+'A'


# In[29]:


#metadata['Phase']=metadata['Phase']+'A'


# In[30]:


metadata


# In[31]:


G = gt.Graph(data_pca, knn=int(7), use_pygsp=True)


# In[32]:


meld_op = meld.MELD(beta=67)
sample_densities = meld_op.fit_transform(G, sample_labels=metadata['sample2'])


# In[33]:


sample_densities


# In[34]:


sample_likelihoods = sample_densities.copy()


# In[35]:


sample_likelihoods


# In[36]:


curr_cols = sample_densities.columns[[col.endswith('A') for col in sample_densities.columns]]
sample_likelihoods[curr_cols] = sklearn.preprocessing.normalize(sample_densities[curr_cols], norm='l1')


# In[37]:


sample_likelihoods


# In[38]:


experimental_samples = ['IR_6hA']
curr_sample = experimental_samples[0]


# In[39]:


curr_sample


# In[40]:


scprep.plot.scatter2d(data_phate, c=sample_likelihoods[curr_sample],
                          vmin=0, vmax=1,
                          title=curr_sample, ticks=False)


# In[41]:


fig, axes = plt.subplots(1,2, figsize=(8.7,4))

scprep.plot.scatter2d(data_phate, c=sample_likelihoods[experimental_samples].mean(axis=1), 
                      cmap=meld.get_meld_cmap(), vmin=0, vmax=1,
                      title='Mean', ticks=False, ax=axes[0])
scprep.plot.scatter2d(data_phate, c=sample_likelihoods[experimental_samples].std(axis=1), vmin=0, 
                      cmap='inferno', title='St. Dev.', ticks=False, ax=axes[1])

fig.tight_layout()


# In[42]:


metadata['chd_likelihood'] = sample_likelihoods[experimental_samples].mean(axis=1).values


# In[49]:


metadata['clusterID'] = scprep.utils.sort_clusters_by_values(metadata['comb.cluster'], metadata['chd_likelihood'])


# In[44]:


#metadata['clusterID'] = scprep.utils.sort_clusters_by_values(metadata['Phase'], metadata['chd_likelihood'])


# In[45]:


#metadata['Phase']= metadata['PhaseA'].replace('A','',regex=True)


# In[46]:


metadata['Phase']


# In[50]:


metadata


# In[52]:


type(metadata)


# In[56]:


metadata.to_csv("/home/liyaru/res/result/6h_metadata.csv",sep=",",index=True,header=True)


# In[47]:


#xaxis group yaxis chd color sample

fig, ax = plt.subplots(1, figsize=(15,10))

# See example usage: https://scprep.readthedocs.io/en/stable/examples/jitter.html
scprep.plot.jitter(metadata['clusterID'], metadata['chd_likelihood'], c=metadata['sample2'], 
                   legend=False, plot_means=False, xlabel=False, ylabel='Mean chd likelihood',
                   ax=ax)

### This code will plot the ratio of tyr:chd cells per cluster
means = metadata.groupby('clusterID')['IR'].mean()
ax.scatter(means.index, means - np.mean(metadata['IR']) + 0.5, edgecolor='k', s=100)

# Axis tick labels
ax.set_xticklabels(metadata.set_index('clusterID')[''].drop_duplicates().sort_index(), rotation=90)
ax.set_ylim(0,1)

fig.tight_layout()


# In[48]:


#xaxis group yaxis chd color sample

fig, ax = plt.subplots(1, figsize=(15,10))

# See example usage: https://scprep.readthedocs.io/en/stable/examples/jitter.html
scprep.plot.jitter(metadata['clusterID'], metadata['chd_likelihood'], c=metadata['sample2'], 
                   legend=False, plot_means=False, xlabel=False, ylabel='Mean chd likelihood',
                   ax=ax)

### This code will plot the ratio of tyr:chd cells per cluster
means = metadata.groupby('clusterID')['IR'].mean()
ax.scatter(means.index, means - np.mean(metadata['IR']) + 0.5, edgecolor='k', s=100)

# Axis tick labels
ax.set_xticklabels(metadata.set_index('clusterID')['comb.cluster'].drop_duplicates().sort_index(), rotation=90)
ax.set_ylim(0,1)

fig.tight_layout()


# In[ ]:


metadata


# In[ ]:


# Get cluster indicides and number of cells per cluster
clusters, counts = np.unique(metadata['clusterID'], return_counts=True)

# Keep cluster labels with at least 1% of the data
clusters = clusters[counts > data.shape[0] * 0.01]


# In[ ]:


data_cluster_phate = {}

for cluster in clusters:
    curr_data = data_pca.loc[metadata['clusterID'] == cluster]
    data_cluster_phate[cluster] = phate.PHATE(verbose=0).fit_transform(curr_data)


# In[ ]:


fig,axes= plt.subplots(5,4, figsize=(4*3, 5*3))

for i , ax in enumerate(axes.flatten()):
    if not i < len(clusters):
        ax.axis('off')
        continue
    curr_cluster = clusters[i]
    curr_phate = data_cluster_phate[curr_cluster]
    
    scprep.plot.scatter2d(curr_phate, 
                          c=metadata['chd_likelihood'].loc[metadata['clusterID'] == curr_cluster], 
                          cmap=meld.get_meld_cmap(), vmin=0, vmax=1,
                         ax=ax, ticks=False, 
                          title='Cluster {} ({})'.format(curr_cluster, curr_phate.shape[0]), 
                          legend=False, fontsize=10)


# In[ ]:


np.random.seed(0)
vfc_op_per_cluster = {}

for cluster in np.unique(metadata['clusterID']):
    curr_G = gt.Graph(data_pca.loc[metadata['clusterID'] == cluster], use_pygsp=True)
    curr_G.compute_fourier_basis()
    curr_sample_labels = metadata['IR'].loc[metadata['clusterID'] == cluster]
    curr_likelihood = metadata['chd_likelihood'].loc[metadata['clusterID'] == cluster]
    curr_vfc = meld.VertexFrequencyCluster(n_clusters = 3)
    curr_vfc.fit_transform(curr_G, curr_sample_labels, curr_likelihood)
    vfc_op_per_cluster[cluster] = curr_vfc


# In[ ]:


subclustering_results = {}
for cluster in np.unique(metadata['clusterID']):
    curr_vfc = vfc_op_per_cluster[cluster]
    clusters_by_n = {}
    for n in [2,3,4,5]:
        clusters_by_n[n] = curr_vfc.predict(n)
    subclustering_results[cluster] = clusters_by_n

