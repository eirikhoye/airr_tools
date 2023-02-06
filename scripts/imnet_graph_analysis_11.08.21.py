#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import scipy
import networkx as nx
import igraph as ig


# In[2]:


get_ipython().system('pwd')


# In[3]:


df_dict = {}
path_dir = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/immunarch_data_376276/'
for file in os.listdir(path_dir):
    if not file.startswith('LivMet') and not file.endswith('.tsv'):
        continue
    dataset = pd.read_csv(path_dir+file, sep='\t')
    sample_name = dataset.sample_name.unique()[0]
    df_dict[sample_name] = dataset
    
path_dir = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/immunarch_data_12511620/'
for file in os.listdir(path_dir):
    if not file.startswith('LivMet') and not file.endswith('.tsv'):
        continue
    dataset = pd.read_csv(path_dir+file, sep='\t')
    sample_name = dataset.sample_name.unique()[0]
    df_dict[sample_name] = dataset

    
df_rearrangement = pd.concat(df_dict.values(), ignore_index=True)
df_rearrangement

public_clones = pd.read_csv('/Users/eirikhoy/Dropbox/projects/airr_comet/data/public_clones_sharing_level.csv')
df_rearrangement = pd.merge(df_rearrangement, public_clones, how='left')
df_rearrangement['sharing_level'] = df_rearrangement['sharing_level'].fillna(0)
df_rearrangement


# In[4]:


path_dir = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/imnet_output/rearrangement_files_with_sharing_all/'
g_n_dict = {}
graph_dict = {}
for file in os.listdir(path_dir):
    if not file.endswith('.graphml'):
        continue
    sample_name = file.split('.')[0]
    g = nx.read_graphml(path_dir+file)
    graph_dict[sample_name] = g
    g_n_dict[sample_name] = {'nr_nodes' : g.number_of_nodes(), 'nr_edges' : g.number_of_edges()}
df = pd.DataFrame.from_dict(g_n_dict).transpose()
df = df.reset_index().rename(columns={'index':'Sample'})
df


# In[5]:


metadata = pd.concat([
    pd.read_csv('/Users/eirikhoy/Dropbox/projects/airr_comet/data/immunarch_data_376276/metadata.txt', sep='\t'), 
    pd.read_csv('/Users/eirikhoy/Dropbox/projects/airr_comet/data/immunarch_data_12511620/metadata.txt', sep='\t')])
metadata


# ### Mean Sharing Level

# In[6]:


df_rearrangement_meta = pd.merge(df_rearrangement, metadata, left_on='sample_name', right_on='Sample')
df_rearrangement_meta

#df_rearrangement_meta['mean_share'] = df_rearrangement_meta.groupby('group').transform('mean')['sharing_level']
#df_rearrangement_meta#sns.barplot(x='group', y='mean_share', data=df_rearrangement_meta[df_rearrangement_meta.group != 'na'])def public_private(x):
    if x == 0:
        return(0)
    if x != 0:
        return(1)
    
df_rearrangement_meta['public_private'] = list(map(public_private, df_rearrangement_meta['sharing_level']))
df_rearrangement_meta['unique_clone'] = 1
df_rearrangement_metadf_rearrangement_meta['total_public'] = df_rearrangement_meta.groupby('Sample')['public_private'].transform('sum')
df_rearrangement_meta['total_clones'] = df_rearrangement_meta.groupby('Sample')['unique_clone'].transform('sum')
df_rearrangement_meta['percentage_clones'] = df_rearrangement_meta.total_public / df_rearrangement_meta.total_clones * 100
df_rearrangement_metasns.barplot(data=df_rearrangement_meta[df_rearrangement_meta.group != 'na'][['sample_name', 'group', 'percentage_clones']].drop_duplicates()
, x='group', y='percentage_clones')
# When looking at public clones defined by aa sequences present in more than one dataset across all datasets, there is no difference between groups

# In[7]:


df_rearrangement_meta#.sharing_level.unique()


# In[8]:


df['Sample_short'] = [d.split('_group')[0] for d in df.Sample]
df_meta = pd.merge(df, metadata, left_on='Sample_short', right_on='COMET_ID')
df_meta['group'] = pd.Categorical(df_meta['group'], ['no_chemo', 'short_chemo', 'long_chemo'])
df_meta = df_meta.sort_values(by='group')
df_meta = df_meta.dropna()
df_meta


# In[9]:


sns.set(rc={'figure.figsize':(8,6)})
sns.set_theme(style='whitegrid', palette='muted')

ax = sns.swarmplot(data=df_meta[df_meta.group != 'na'], x='group', y='nr_nodes', hue='Kit_ID')


# Slight increase in short chemo for number of unique clones

# In[10]:


sns.set(rc={'figure.figsize':(12,4)})
sns.set_theme(style='whitegrid', palette='muted')

ax = sns.boxplot(data=df_meta, y='group', x='nr_edges', whis=np.inf, palette='Dark2')
ax = sns.swarmplot(data=df_meta, y='group', x='nr_edges', color='.2')


# In[11]:


sns.set(rc={'figure.figsize':(8,6)})
sns.set_theme(style='whitegrid', palette='muted')

ax = sns.barplot(data=df_meta[df_meta.group != 'na'], x='group', y='nr_edges',
                palette='viridis')


# Also slightly higher number of edges in the short chemo group

# In[12]:


sns.set_style('white')
ax = sns.scatterplot(x='nr_nodes', y='nr_edges', data=df_meta[df_meta.group != 'na'], color='black')
sns.despine()
ax.set(xlabel="Number of Clonotypes (Nodes)", ylabel='Number of Connections (Edges)')


# In[13]:


df_meta


# Looks like connectivity is clearly a function of the number of unique clones

# ### Only Kit 1

# In[14]:


sns.set(rc={'figure.figsize':(8,6)})
sns.set_theme(style='whitegrid', palette='muted')

ax = sns.swarmplot(data=df_meta[(df_meta.group != 'na') & (df_meta.Kit_ID == 376276)], x='group', y='nr_nodes', hue='Kit_ID')


# In[ ]:





# In[15]:


sns.set(rc={'figure.figsize':(8,6)})
sns.set_theme(style='whitegrid', palette='muted')

ax = sns.swarmplot(data=df_meta[(df_meta.group != 'na') & (df_meta.Kit_ID == 376276)], x='group', y='nr_edges', hue='Kit_ID')


# In[16]:


sns.scatterplot(x='nr_nodes', y='nr_edges', data=df_meta[(df_meta.group != 'na') & (df_meta.Kit_ID == 376276)])


# In[17]:


sns.set(rc={'figure.figsize':(8,6)})
sns.set_theme(style='whitegrid', palette='muted')

ax = sns.swarmplot(data=df_meta[(df_meta.group != 'na') & (df_meta.Kit_ID == 12511620)], x='group', y='nr_nodes', hue='Kit_ID')


# In[18]:


df_meta.sort_values(by='nr_nodes', ascending=False)


# In[19]:


sns.set(rc={'figure.figsize':(8,6)})
sns.set_theme(style='whitegrid', palette='muted')

ax = sns.swarmplot(data=df_meta[(df_meta.group != 'na') & (df_meta.Kit_ID == 12511620)], x='group', y='nr_edges', hue='Kit_ID')


# In[20]:


sns.scatterplot(x='nr_nodes', y='nr_edges', data=df_meta[(df_meta.group != 'na') & (df_meta.Kit_ID == 12511620)])


# In[21]:


df_meta.sort_values(by='group')


# In[22]:


df_meta.nr_nodes.mean()


# In[23]:


df_meta.nr_nodes.max()


# In[24]:


df


# In[25]:


max_comp_dict = {}
for sample_name in graph_dict.keys():
    max_comp_dict[sample_name] = [np.array([len(c) for c in nx.connected_components(graph_dict[sample_name])]).max(),
                                  len(graph_dict[sample_name]),
                                 np.array([len(c) for c in nx.connected_components(graph_dict[sample_name])]).max() / len(graph_dict[sample_name]) * 100]
                                  


# In[26]:


max_comp_df = pd.DataFrame.from_dict(max_comp_dict).transpose()
max_comp_df.columns = ['max_clust', 'nr_sequences', 'percentage']
max_comp_df = max_comp_df.reset_index()#.rename(columns={'Samples':'index'})
max_comp_df


# In[27]:


pd.merge(df, metadata, on='Sample')


# In[28]:


sns.set(rc={'figure.figsize':(4,6)})
sns.set_theme(style='whitegrid', palette='muted')

ax = sns.swarmplot(data=max_comp_df, y='percentage', color='.2')


# In[29]:


sns.scatterplot(data=max_comp_df, x='nr_sequences', y='percentage')


# In[30]:


max_comp_df.corr()


# In[31]:


ig


# In[32]:


path_dir = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/imnet_output/graphs/'
g = ig.Graph.Read_GraphML(path_dir+'LivMet_1.tsv.graphml')


# In[33]:


path_dir = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/imnet_output/graphs/'
graph_dict = {}
for file in os.listdir(path_dir):
    if not file.endswith('tsv.graphml'):
        continue
    sample_name = file.split('.')[0]
    g = ig.Graph.Read_GraphML(path_dir+file)
    graph_dict[sample_name] = g
graph_dict


# In[34]:


kit_1_qc = pd.read_csv('/Users/eirikhoy/Dropbox/projects/airr_comet/data/qc_reports/qcReport_kit_1.tsv', sep='\t', index_col=False)
kit_2_qc = pd.read_csv('/Users/eirikhoy/Dropbox/projects/airr_comet/data/qc_reports/qcReport_kit_2.tsv', sep='\t', index_col=False)
qc_report = pd.concat([kit_1_qc, kit_2_qc])
qc_report = qc_report[~qc_report['Sample Name'].isin(['Positive', 'Negative'])]
qc_report


# In[35]:


low_coverage_samples = list(qc_report[qc_report.Coverage < 5]['Sample Name'])
low_coverage_samples


# In[36]:


is_first=True
for dataset in graph_dict:
    
    global_param_ = {
        'nr_nodes'          : len(graph_dict[dataset].vs),
        'number_edges'      : len(graph_dict[dataset].es),
        'largest_component' : max([len(i) for i in graph_dict[dataset].clusters()]),
        'maximal_k_core'    : max(graph_dict[dataset].coreness()),
        'diameter'          : graph_dict[dataset].diameter(),
        'assortativity'     : graph_dict[dataset].assortativity_degree(directed=False),
        'average_cluster'   : np.mean([len(i) for i in graph_dict[dataset].clusters()]),
        'average_degree'    : np.mean(sorted(graph_dict[dataset].indegree(), reverse=True)),
        'clustering_coeff'  : np.mean(graph_dict[dataset].transitivity_local_undirected(mode='zero')),
        'density'           : np.mean(graph_dict[dataset].density(loops=False))
    }
    if is_first:
        global_param = pd.DataFrame.from_dict(global_param_, orient='index', columns=[dataset])
        is_first=False
    else:
        global_param = pd.concat([global_param, pd.DataFrame.from_dict(global_param_, orient='index', columns=[dataset])],
                                axis=1)
global_param = global_param.transpose()
global_param_subset = global_param.drop(low_coverage_samples, axis=0)

global_param_subset = global_param_subset.drop(['LivMet_38', 'LivMet_11', 'LivMet_17',
                   'LivMet_40', 'LivMet_42', 'LivMet_33', 'LivMet_21'], axis=0)
global_param_subset


# In[37]:


global_param_subset.to_csv(path_or_buf='/Users/eirikhoy/Dropbox/projects/airr_comet/data/global_param_subset.tsv', sep='\t')


# In[38]:


global_param_subset.describe()


# In[39]:


import scipy.stats


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h


# In[40]:


mean_confidence_interval(global_param_subset.nr_nodes)


# In[41]:


mean_confidence_interval(global_param_subset.number_edges)


# In[42]:


global_param_subset['fraction'] = global_param_subset.number_edges / global_param_subset.nr_nodes


# In[43]:


mean_confidence_interval(global_param_subset.fraction)


# In[44]:


sns.violinplot(y='nr_nodes', data=global_param_subset)
sns.swarmplot(y='nr_nodes', data=global_param_subset, size=5, color='black')


# In[45]:


is_first=True
for dataset in graph_dict:

    local_param_ = {
        'transitivity'           : np.mean(graph_dict[dataset].transitivity_local_undirected(mode='zero')),
        'authority'              : np.mean(graph_dict[dataset].authority_score()),
        'pagerank'               : np.mean(graph_dict[dataset].personalized_pagerank(directed=False)),
        'eigenvector_centrality' : np.mean(graph_dict[dataset].eigenvector_centrality(directed=False)),
        'closeness'              : np.mean(graph_dict[dataset].closeness()),
        'betweenness'            : np.mean(graph_dict[dataset].betweenness(directed=False))
    }
    if is_first:
        local_param = pd.DataFrame.from_dict(local_param_, orient='index', columns=[dataset])
        is_first=False
    else:
        local_param = pd.concat([local_param, pd.DataFrame.from_dict(local_param_, orient='index', columns=[dataset])],
                               axis=1)

local_param = local_param.transpose()
#local_param = local_param.reset_index().rename(columns={'index':'Sample'})
#local_param = pd.merge(local_param, metadata, on='Sample')
local_param


# In[46]:


local_param.to_csv(path_or_buf='/Users/eirikhoy/Dropbox/projects/airr_comet/data/local_param.tsv', sep='\t')


# In[47]:


local_param.describe()


# In[48]:


# Get list of networks with few nodes
low_nodes = list(global_param[global_param.nr_nodes < 5000].index)
low_nodes


# ### Global Parameters

# In[49]:


sns.set(rc={'figure.figsize':(30,4)})
sns.set_theme(style='whitegrid', palette='muted')
global_param_norm = global_param
#lobal_param_norm.index = global_param_norm.Sample
#global_param_norm = global_param_norm.drop('Sample', axis=1)
global_param_norm = (global_param_norm-global_param_norm.mean()) / global_param_norm.std()
sns.heatmap(global_param_norm.transpose(), cmap='viridis')
global_param_norm


# #### Global Parameters drop small datasets

# In[50]:


pd.merge(global_param_norm, metadata[['Sample', 'group']], left_on=global_param_norm.index, right_on='Sample')


# In[51]:


sns.set(rc={'figure.figsize':(30,4)})
sns.set_theme(style='whitegrid', palette='muted')
global_param_norm = global_param.drop(low_nodes, axis=0)
#global_param_norm = pd.merge(global_param_norm, metadata[['Sample', 'group']], left_on=global_param_norm.index, right_on='Sample')
#lobal_param_norm.index = global_param_norm.Sample
#global_param_norm = global_param_norm.drop('Sample', axis=1)
global_param_norm = (global_param_norm-global_param_norm.mean()) / global_param_norm.std()
sns.heatmap(global_param_norm.transpose(), cmap='viridis')
global_param_norm


# In[ ]:





# In[52]:


sns.set(rc={'figure.figsize':(30,4)})
sns.set_theme(style='whitegrid', palette='muted')
global_param_norm = global_param.drop(low_nodes, axis=0)
#global_param_norm = pd.merge(global_param_norm, metadata[['Sample', 'group']], left_on=global_param_norm.index, right_on='Sample')
#lobal_param_norm.index = global_param_norm.Sample
#global_param_norm = global_param_norm.drop('Sample', axis=1)
global_param_norm = (global_param_norm-global_param_norm.mean()) / global_param_norm.std()
sns.clustermap(global_param_norm, cmap='viridis')
global_param_norm


# In[53]:


sns.set(rc={'figure.figsize':(10,8)})

sns.heatmap(global_param.drop(low_nodes, axis=0).corr(), vmin=0, vmax=1)


# In[54]:


global_param


# ### Local Parameters

# In[55]:


local_param.describe()


# In[56]:


sns.set(rc={'figure.figsize':(30,4)})
sns.set_theme(style='whitegrid', palette='muted')
local_param_norm = local_param
local_param_norm.index = local_param_norm.Sample
local_param_norm = local_param_norm.drop('Sample', axis=1)
local_param_norm = (local_param_norm-local_param_norm.mean()) / local_param_norm.std()
sns.heatmap(local_param_norm.transpose(), cmap='viridis')
local_param_norm


# #### Local Parameters drop small datasets

# In[57]:


sns.set(rc={'figure.figsize':(30,4)})
sns.set_theme(style='whitegrid', palette='muted')
local_param_norm = local_param.drop(low_nodes, axis=0)
local_param_norm.index = local_param_norm.Sample
local_param_norm = local_param_norm.drop('Sample', axis=1)
local_param_norm = (local_param_norm-local_param_norm.mean()) / local_param_norm.std()
sns.heatmap(local_param_norm.transpose(), cmap='viridis')
local_param_norm


# ### Comparison Groups Global

# In[58]:


#global_param_meta =
global_param_meta = global_param.copy()
global_param_meta['Sample'] = global_param_meta.index
global_param_meta = global_param_meta.reset_index(drop=True)
global_param_meta = pd.merge(global_param_meta, metadata, how='inner')
global_compare = global_param_meta[global_param_meta.group != 'na'].groupby('group').mean().drop(['ID', 'COMET_ID_nr', 'Kit_ID'], axis=1)
global_compare


# In[59]:


global_compare_norm = (global_compare - global_compare.mean()) / global_compare.std()
global_compare_norm


# In[60]:


sns.set(rc={'figure.figsize':(10,4)})
sns.heatmap(global_compare_norm, cmap='viridis')


# In[61]:


global_param_meta['group'] = pd.Categorical(global_param_meta['group'], ['no_chemo', 'short_chemo', 'long_chemo'])
global_param_meta = global_param_meta.sort_values('group')
global_param_meta


# In[62]:


sns.set_context(rc={'patch.linewidth':1.0})
plt.figure(figsize=(4, 8))

fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(16, 8))

sns.barplot(x='group', y='nr_nodes', data=global_param_meta[global_param_meta.group != 'na'], palette='GnBu', edgecolor='black', ax=ax1)
sns.barplot(x='group', y='number_edges', data=global_param_meta[global_param_meta.group != 'na'], palette='GnBu', edgecolor='black', ax=ax2)
sns.barplot(x='group', y='average_degree', data=global_param_meta[global_param_meta.group != 'na'], palette='GnBu', edgecolor='black', ax=ax3)
sns.barplot(x='group', y='largest_component', data=global_param_meta[global_param_meta.group != 'na'], palette='GnBu', edgecolor='black', ax=ax4)
fig.tight_layout()


# In[63]:


path_dir = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/imnet_output_sharing_level/graphs/'
graph_dict = {}
for file in os.listdir(path_dir):
    if not file.endswith('.graphml'):
        continue
    sample_name = file.split('.')[0]
    g = ig.Graph.Read_GraphML(path_dir+file)
    graph_dict[sample_name] = g
graph_dict


# In[64]:


len(graph_dict['LivMet_24'].indegree())


# In[65]:


is_first = True
for dataset in graph_dict:
    if is_first:
        degree_share = pd.DataFrame({
            'degree' : graph_dict[dataset].indegree(),
            'share_level' : graph_dict[dataset].vs['share_level']
        })
        is_first = False
    else:
        degree_share_ = pd.DataFrame({
            'degree' : graph_dict[dataset].indegree(),
            'share_level' : graph_dict[dataset].vs['share_level']
        })
        degree_share = pd.concat([degree_share, degree_share_])
        
degree_share


# In[66]:


degree_share.to_csv(path_or_buf='/Users/eirikhoy/Dropbox/projects/airr_comet/data/degree_share.tsv', sep='\t', index=False)


# In[67]:


sns.set_style('white')
plt.figure(figsize=(16, 12))

ax = sns.boxplot(x='degree', y='share_level', data=degree_share, color='lightgrey', fliersize=0.5)
sns.despine()
ax.set(xlabel='Clonal connections (Degrees)', ylabel='Number of patients')


# In[68]:


X = degree_share.sort_values(by='degree').degree.values
y = degree_share.sort_values(by='degree').share_level.values
X = X.reshape(-1,1)
y = y.reshape(-1,1)


# In[69]:


from sklearn.linear_model import LinearRegression
lin_reg = LinearRegression()
lin_reg.fit(X, y)
lin_reg.intercept_, lin_reg.coef_


# In[70]:


max(X)


# In[71]:


y_pred = lin_reg.predict(X)


# In[72]:


len(y_pred)


# In[73]:


from sklearn.metrics import mean_squared_error
mean_squared_error(y, y_pred)


# In[74]:


sns.set_style('white')
plt.figure(figsize=(12, 7))

ax = sns.boxplot(x='degree', y='share_level', data=degree_share, color='lightgrey', fliersize=0.5)
sns.despine()
ax.set(xlabel='Clonal connections (Degrees)', ylabel='Number of patients')
plt.plot(X, y_pred, 'r-', linewidth=2, label='patients = {:.2f} + {:.2f} * degrees'.format(lin_reg.intercept_[0], lin_reg.coef_[0, 0]))
plt.legend(loc='upper left', fontsize=14)
plt.xlim(0, 48)


# In[75]:


degree_share = degree_share.astype('int32')


# In[1]:


degree_share


# In[76]:


sns.set_style('white')
plt.figure(figsize=(9, 7))

ax = sns.barplot(x='share_level', y='degree', data=degree_share, color='black')
#ax.set(ylabel='Clonal connections (Degrees)', xlabel='Number of patients')
plt.ylabel('Clonal connections (Degrees)', fontsize=16)
plt.xlabel('Number of patients', fontsize=16)

ax.set_xticks([0, 10, 20, 30, 40])
ax.set_xticklabels([0, 10, 20, 30, 40])
sns.despine()


# In[77]:


ax.get_xticks()


# In[78]:


from sklearn import metrics
print('Mean Absolute Error:', metrics.mean_absolute_error(y, y_pred))  
print('Mean Squared Error:', metrics.mean_squared_error(y, y_pred))  
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y, y_pred)))


# In[79]:


metrics.r2_score(y, y_pred)


# In[80]:


import statsmodels.api as sm
from scipy import stats

X2 = sm.add_constant(X)
X2
est = sm.OLS(y, X2)
est2 = est.fit()
print(est2.summary2())







# In[81]:


est2.params


# In[82]:


degree_clone


# In[83]:


X = degree_clone.sort_values(by='degree').degree.values
y = degree_clone.sort_values(by='degree').clone_size.values
X = X.reshape(-1,1)
y = y.reshape(-1,1)


# In[84]:


from sklearn.linear_model import LinearRegression
lin_reg = LinearRegression()
lin_reg.fit(X, y)
lin_reg.intercept_, lin_reg.coef_


# In[85]:


y_pred = lin_reg.predict(X)


# In[86]:


sns.set_style('white')
ax = sns.boxplot(x='degree', y='clone_size', data=degree_clone, palette='viridis')
sns.despine()
ax.set(xlabel='Clonal connections (Degrees)', ylabel='Number of T cells in Clonotype')
plt.plot(X, y_pred, 'r-', linewidth=2, label='clonal size = {:.2e} + {:.2e} * degrees'.format(lin_reg.intercept_[0], lin_reg.coef_[0, 0]))
plt.legend(loc='upper left', fontsize=14)
#plt.xlim(0, 48)


# In[87]:


sns.set_style('white')
ax = sns.scatterplot(x='nr_nodes', y='nr_edges', data=df_meta[df_meta.group != 'na'], color='black')
sns.despine()
ax.set(xlabel="Number of Clonotypes", ylabel='Number of Edges')



# In[88]:


# alternatively f_regression from sklearn.feature_selection could be used
# Since the p-values are obtained through certain statistics, we need the 'stat' module from scipy.stats
import numpy as np
import pandas as pd
from sklearn import linear_model
import scipy.stats as stat

# Since we are using an object oriented language such as Python, we can simply define our own 
# LinearRegression class (the same one from sklearn)
# By typing the code below we will ovewrite a part of the class with one that includes p-values
# Here's the full source code of the ORIGINAL class: https://github.com/scikit-learn/scikit-learn/blob/7b136e9/sklearn/linear_model/base.py#L362


class LinearRegression(linear_model.LinearRegression):
    """
    LinearRegression class after sklearn's, but calculate t-statistics
    and p-values for model coefficients (betas).
    Additional attributes available after .fit()
    are `t` and `p` which are of the shape (y.shape[1], X.shape[1])
    which is (n_features, n_coefs)
    This class sets the intercept to 0 by default, since usually we include it
    in X.
    """
    
    # nothing changes in __init__
    def __init__(self, fit_intercept=True, normalize=False, copy_X=True,
                 n_jobs=1):
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.copy_X = copy_X
        self.n_jobs = n_jobs

    
    def fit(self, X, y, n_jobs=1):
        self = super(LinearRegression, self).fit(X, y, n_jobs)
        
        # Calculate SSE (sum of squared errors)
        # and SE (standard error)
        sse = np.sum((self.predict(X) - y) ** 2, axis=0) / float(X.shape[0] - X.shape[1])
        se = np.array([np.sqrt(np.diagonal(sse * np.linalg.inv(np.dot(X.T, X))))])

        # compute the t-statistic for each feature
        self.t = self.coef_ / se
        # find the p-value for each feature
        self.p = np.squeeze(2 * (1 - stat.t.cdf(np.abs(self.t), y.shape[0] - X.shape[1])))
        return self


# In[89]:


df_meta['edge_node_fraction'] = df_meta['nr_edges'] / df_meta['nr_nodes']
df_meta


# In[90]:


df_meta.to_csv(path_or_buf='/Users/eirikhoy/Dropbox/projects/airr_comet/')


# In[91]:


nr_nodes = df_meta.sort_values(by='nr_nodes').nr_nodes.values
nr_edges = df_meta.sort_values(by='nr_nodes').nr_edges.values
X = nr_nodes.reshape(-1,1)#6 * np.random.rand(m, 1) - 3
y = nr_edges.reshape(-1,1)#0.5 * X**2 + X + 2 + np.random.randn(m, 1)

from sklearn.preprocessing import PolynomialFeatures
poly_features = PolynomialFeatures(degree=2, include_bias=False)
X_poly = poly_features.fit_transform(X)

lin_reg = LinearRegression()
lin_reg.fit(X_poly,  y)
lin_reg.intercept_, lin_reg.coef_

#x_mod = np.linspace(min(X), max(X), len(X)).reshape(len(X),1)
x_mod_poly = poly_features.transform(X)
y_mod = lin_reg.predict(x_mod_poly)
plt.figure(figsize=(12, 7))

plt.plot(X, y , 'k.')
plt.plot(X, y_mod, 'r-', linewidth=2, label='y={:.2f}+{:.2e}x^2+{:.2e}x'.format(lin_reg.intercept_[0], lin_reg.coef_[0, 0], lin_reg.coef_[0,1]))
plt.xlabel("Number of clonotypes (nodes)")
plt.ylabel("Number of connections (edges)")
plt.legend(loc='upper left', fontsize=14)
#plt.axis([min(X-5000), max(X+5000), -5000, 45000])
plt.show()
print(lin_reg.intercept_[0], lin_reg.coef_[0,0], lin_reg.coef_[0,1])


# In[92]:


lin_reg.p


# In[107]:


nr_nodes = df_meta.sort_values(by='nr_nodes').nr_nodes.values
fraction = df_meta.sort_values(by='nr_nodes').edge_node_fraction.values 
X = nr_nodes.reshape(-1,1)#6 * np.random.rand(m, 1) - 3
y = fraction.reshape(-1,1)#0.5 * X**2 + X + 2 + np.random.randn(m, 1)

from sklearn.preprocessing import PolynomialFeatures
poly_features = PolynomialFeatures(degree=2, include_bias=False)
X_poly = poly_features.fit_transform(X)

lin_reg = LinearRegression()
lin_reg.fit(X,  y)
lin_reg.intercept_, lin_reg.coef_

#x_mod = np.linspace(min(X), max(X), len(X)).reshape(len(X),1)
#x_mod_poly = poly_features.transform(X)
y_mod = lin_reg.predict(X)
plt.figure(figsize=(12, 7))
sns.set_style('white')


plt.plot(X, y , 'k.')
plt.plot(X, y_mod, 'r-', linewidth=2, label='y={:.2f}+{:.2e}x'.format(lin_reg.intercept_[0], lin_reg.coef_[0,0]))
plt.xlabel("Number of clonotypes (nodes)")
plt.ylabel("Connectivity fraction")
plt.legend(loc='upper left', fontsize=14)
sns.despine()

#plt.axis([min(X-5000), max(X+5000), -5000, 45000])
plt.show()
print(lin_reg.intercept_[0], lin_reg.coef_[0,0])
print(metrics.r2_score(y, y_mod))


# In[113]:


pd.DataFrame({
'nr_nodes' : df_meta.sort_values(by='nr_nodes').nr_nodes.values,
'fraction' : df_meta.sort_values(by='nr_nodes').edge_node_fraction.values 
}).to_csv('../data/connectivity_fraction_to_nr_nodes.tsv', sep='\t')


# In[94]:


from sklearn import metrics

print('Mean Absolute Error:', metrics.mean_absolute_error(y, y_pred))  
print('Mean Squared Error:', metrics.mean_squared_error(y, y_pred))  
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y, y_pred)))


# In[95]:


a = lin_reg.intercept_
b = lin_reg.coef_[0, 0]
c = lin_reg.coef_[0, 1]
popt = np.array([a, b, c])


# In[ ]:





# In[96]:


popt


# In[97]:


df_meta.head()


# In[98]:


path_dir = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/evenness_auc.tsv'
evenness_auc = pd.read_csv(path_dir, sep='\t')
evenness_auc.head()


# In[99]:


evenness_auc


# In[ ]:





# In[100]:


df_meta


# In[101]:


edges_even_auc = pd.merge(df_meta[['Sample_y', 'nr_edges','nr_nodes']], evenness_auc[['Sample', 'Evenenss_AUC']],
         left_on='Sample_y', right_on='Sample').drop('Sample_y', axis=1)
edges_even_auc['fraction'] = edges_even_auc.nr_edges / edges_even_auc.nr_nodes
edges_even_auc.head()


# In[102]:


sns.set_style('white')
plt.figure(figsize=(8, 6))
ax = sns.scatterplot(x='Evenenss_AUC', y='fraction', data=edges_even_auc, color='black')
sns.despine()

ax.set(xlabel="Monoclonality (Evenness AUC)", ylabel='Edges to nodes fraction')


# In[114]:


X = edges_even_auc.sort_values(by='Evenenss_AUC').Evenenss_AUC.values
y = edges_even_auc.sort_values(by='Evenenss_AUC').fraction.values
X = X.reshape(-1,1)
y = y.reshape(-1,1)

#from sklearn.linear_model import LinearRegression
lin_reg = LinearRegression()
lin_reg.fit(X, y)
print(lin_reg.intercept_, lin_reg.coef_)

y_pred = lin_reg.predict(X)

from sklearn.metrics import mean_squared_error
print(mean_squared_error(y, y_pred))

sns.set_style('white')
plt.figure(figsize=(12, 7))
ax = sns.scatterplot(x='Evenenss_AUC', y='fraction', data=edges_even_auc, color='black')
sns.despine()

ax.set(xlabel="Monoclonality (Evenness AUC)", ylabel='Connectivity fraction')
plt.plot(X, y_pred, 'r-', linewidth=2, label='patients = {:.2f} + {:.2f} * degrees'.format(lin_reg.intercept_[0], lin_reg.coef_[0, 0]))
plt.legend(loc='upper left', fontsize=14)
#plt.xlim(0, 48)
print(lin_reg.intercept_[0], lin_reg.coef_[0,0])
print(metrics.r2_score(y, y_pred))


# In[116]:


pd.DataFrame({
'evenness_auc' : edges_even_auc.sort_values(by='Evenenss_AUC').Evenenss_AUC.values,
'connectivity_fraction' : edges_even_auc.sort_values(by='Evenenss_AUC').fraction.values
}).to_csv("../data/connectivity_fraction_to_evenness_auc.tsv", sep="\t")


# In[105]:


import statsmodels.api as sm
from scipy import stats

X2 = sm.add_constant(X)
X2
est = sm.OLS(y, X2)
est2 = est.fit()
print(est2.summary2())







# In[96]:


edge_node = df_meta[['Sample_y', 'nr_edges', 'nr_nodes']]
edge_node['edge_to_node'] = edge_node.nr_edges / edge_node.nr_nodes
edge_node.head()


# In[97]:


edges_even_auc = pd.merge(edge_node[['Sample_y', 'edge_to_node']], evenness_auc[['Sample', 'Evenenss_AUC']],
         left_on='Sample_y', right_on='Sample').drop('Sample_y', axis=1)
edges_even_auc.head()


# In[98]:


sns.scatterplot(x='Evenenss_AUC', y='edge_to_node', data=edges_even_auc)


# In[99]:


sns.scatterplot(x='nr_nodes', y='nr_edges', data=edge_node, color='black')


# In[100]:


sns.histplot(edge_node['edge_to_node'])


# # |
# # |
# # |
# # |
# # |
# # |
# # |
# # |
# # |

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ### Comparison Groups Local

# In[101]:


local_param_meta = pd.merge(local_param.reset_index(drop=True), metadata, how='inner')


# In[102]:



local_compare = local_param_meta[local_param_meta.group != 'na'].groupby('group').mean().drop(['ID', 'COMET_ID_nr', 'Kit_ID'], axis=1)
local_compare_norm = (local_compare - local_compare.mean()) / local_compare.std()
local_compare_norm


# In[103]:


local_param_meta['group'] = pd.Categorical(local_param_meta['group'], ['no_chemo', 'short_chemo', 'long_chemo'])
local_param_meta = local_param_meta.sort_values('group')


# In[104]:


local_param_meta


sns.set_context(rc={'patch.linewidth':1.0})
plt.figure(figsize=(4, 8))

fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(nrows=1, ncols=5, figsize=(20, 8))
sns.barplot(x='group', y='authority', data=local_param_meta[local_param_meta.group != 'na'], palette='GnBu', edgecolor='black', ax=ax1)
sns.barplot(x='group', y='pagerank', data=local_param_meta[local_param_meta.group != 'na'], palette='GnBu', edgecolor='black', ax=ax2)
sns.barplot(x='group', y='eigenvector_centrality', data=local_param_meta[local_param_meta.group != 'na'], palette='GnBu', edgecolor='black', ax=ax3)
sns.barplot(x='group', y='closeness', data=local_param_meta[local_param_meta.group != 'na'], palette='GnBu', edgecolor='black', ax=ax4)
sns.barplot(x='group', y='betweenness', data=local_param_meta[local_param_meta.group != 'na'], palette='GnBu', edgecolor='black', ax=ax5)
fig.tight_layout()


# In[105]:


sns.set(rc={'figure.figsize':(10,4)})
sns.heatmap(local_compare_norm, cmap='viridis')


# In[106]:


df_rearrangement[df_rearrangement.sample_name == 'LivMet_66'].sort_values(by='templates')


# In[83]:


# subsampling 100 times, find s

df_66 = df_rearrangement[df_rearrangement.sample_name == 'LivMet_66']

df_66.sort_values(by='templates')


# In[84]:


test = df_66.sample(2000)
test['frequency'] = test.templates / test.templates.sum()
test.sort_values(by='templates')


# In[ ]:


unsampled_max_clone_size = df_66.frequency.max()
unsampled_max_clone_size


# In[ ]:


#def max_clone_size(df_rearrangement)

subsample_dict = {}
max_clone_size_dict = {}

for i in range(1,101):
    subsample = df_66.sample(2000)
    subsample['frequency'] = subsample.templates / subsample.templates.sum()
    subsample_dict[i] = subsample
    max_clone_size = subsample.frequency.max()
    max_clone_size_dict[i] = max_clone_size
median_subsample_max = np.median(list(max_clone_size_dict.values()))

median_subsample_diff = {}
for i in max_clone_size_dict:
    diff = max_clone_size_dict[i] - median_subsample_max
    median_subsample_diff[i] = diff

best_subsample = subsample_dict[min(median_subsample_diff, key=median_subsample_diff.get)]

best_subsample.sort_values(by='templates')


# In[ ]:


import sys
sys.path.append('../src/')
import graph
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '')

g = graph.RepertoireGraph(best_subsample)
g.levenstein(1)
g.threshold(1)
g.make_graph(1)
p = g.plot_graph(1)


# In[ ]:


p['network']


# In[ ]:


out_path = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/graph_clone_sampling/graphs/'
#ig.write_graph(g._graph[1], file=out_path+'_'+'LivMet_66.graphml', format='graphml')
g._graph[1].write_graphml(out_path+'clone_sample_'+'LivMet_66.graphml')


# In[ ]:


top_5000_clones = df_66.sort_values(by='templates', ascending=False).head(5000)
top_5000_clones


# In[ ]:


import sys
sys.path.append('../src/')
import graph
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '')

g = graph.RepertoireGraph(top_5000_clones)
g.levenstein(1)
g.threshold(1)
g.make_graph(1)
p = g.plot_graph(1)


# In[ ]:


p['network']


# In[ ]:


out_path = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/graph_top_5000/graphs/'
#ig.write_graph(g._graph[1], file=out_path+'_'+'LivMet_66.graphml', format='graphml')
g._graph[1].write_graphml(out_path+'top_5000_'+'LivMet_66.graphml')


# In[ ]:


g.assortativity[1]


# LivMet_8
# 

# In[ ]:


dataset = 'LivMet_8'

df_sample = df_rearrangement[df_rearrangement.sample_name == dataset]
top_5000_clones = df_sample.sort_values(by='templates', ascending=False).head(5000)

import sys
sys.path.append('../src/')
import graph
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '')

g = graph.RepertoireGraph(top_5000_clones)
g.levenstein(1)
g.threshold(1)
g.make_graph(1)
#p = g.plot_graph(1)

out_path = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/graph_top_5000/graphs/'
#ig.write_graph(g._graph[1], file=out_path+'_'+'LivMet_66.graphml', format='graphml')
g._graph[1].write_graphml(out_path+'top_5000_'+dataset+'.graphml')


# LivMet_40
# LivMet_43
# LivMet_44

# In[ ]:


dataset = 'LivMet_40'

df_sample = df_rearrangement[df_rearrangement.sample_name == dataset]
top_5000_clones = df_sample.sort_values(by='templates', ascending=False).head(5000)

import sys
sys.path.append('../src/')
import graph
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '')

g = graph.RepertoireGraph(top_5000_clones)
g.levenstein(1)
g.threshold(1)
g.make_graph(1)
#p = g.plot_graph(1)

out_path = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/graph_top_5000/graphs/'
#ig.write_graph(g._graph[1], file=out_path+'_'+'LivMet_66.graphml', format='graphml')
g._graph[1].write_graphml(out_path+'top_5000_'+dataset+'.graphml')


# In[ ]:


g.n_edges[1]


# In[ ]:


dataset = 'LivMet_44'

df_sample = df_rearrangement[df_rearrangement.sample_name == dataset]
top_5000_clones = df_sample.sort_values(by='templates', ascending=False).head(5000)

import sys
sys.path.append('../src/')
import graph
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '')

g = graph.RepertoireGraph(top_5000_clones)
g.levenstein(1)
g.threshold(1)
g.make_graph(1)
#p = g.plot_graph(1)

out_path = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/graph_top_5000/graphs/'
#ig.write_graph(g._graph[1], file=out_path+'_'+'LivMet_66.graphml', format='graphml')
g._graph[1].write_graphml(out_path+'top_5000_'+dataset+'.graphml')


# In[ ]:


g.n_edges[1]


# In[ ]:


g._graph[1].vs['size']


# In[ ]:


is_first=True
for dataset in graph_dict:

    local_param_ = {
        'transitivity'           : np.mean(graph_dict[dataset].transitivity_local_undirected(mode='zero')),
        'authority'              : np.mean(graph_dict[dataset].authority_score()),
        'pagerank'               : np.mean(graph_dict[dataset].personalized_pagerank(directed=False)),
        'eigenvector_centrality' : np.mean(graph_dict[dataset].eigenvector_centrality(directed=False)),
        'closeness'              : np.mean(graph_dict[dataset].closeness()),
        'betweenness'            : np.mean(graph_dict[dataset].betweenness(directed=False))
    }
    if is_first:
        local_param = pd.DataFrame.from_dict(local_param_, orient='index', columns=[dataset])
        is_first=False
    else:
        local_param = pd.concat([local_param, pd.DataFrame.from_dict(local_param_, orient='index', columns=[dataset])],
                               axis=1)

local_param = local_param.transpose()
#local_param = local_param.reset_index().rename(columns={'index':'Sample'})
#local_param = pd.merge(local_param, metadata, on='Sample')
local_param


# In[ ]:


sharing_level


# In[ ]:


graph_dict['LivMet_37'].vs


# In[ ]:


df_rearrangement[df_rearrangement.sample_name == dataset]


# In[ ]:


outpath = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/imnet_output_sharing_level/graphs/'
for dataset in graph_dict:
    print(dataset)
    clone_share_dict = {k:v for k, v in df_rearrangement[df_rearrangement.sample_name == dataset][['amino_acid', 'sharing_level']].values}
    clone_size_dict = {k:v for k, v in df_rearrangement[df_rearrangement.sample_name == dataset][['amino_acid', 'templates']].values}
    graph_dict[dataset].vs['share_level'] = [clone_share_dict[v] for v in graph_dict[dataset].vs['id']]
    graph_dict[dataset].vs['clone_size']  = [clone_size_dict[v] for v in graph_dict[dataset].vs['id']]
    graph_dict[dataset].write_graphml(outpath+str(dataset)+'.graphml')


# In[ ]:





# In[ ]:





# In[ ]:


[v for v in graph_dict['LivMet_37'].vs]


# In[ ]:





# In[99]:


df_rearrangement_meta


# In[111]:


total = df_rearrangement_meta['templates'].sum()
private = df_rearrangement_meta[df_rearrangement_meta.public_private==0]['templates'].sum()
public = df_rearrangement_meta[df_rearrangement_meta.public_private==1]['templates'].sum()

public_percentage = public / total
public_percentage


# In[113]:


total   = len(df_rearrangement_meta)
private = len(df_rearrangement_meta[df_rearrangement_meta.public_private==0])
public  = len(df_rearrangement_meta[df_rearrangement_meta.public_private==1])

public_percentage = public / total
public_percentage


# In[112]:


len(df_rearrangement_meta)


# In[114]:


df_rearrangement_meta[df_rearrangement_meta['sample_name'] == 'LivMet_35']


# In[ ]:




