# random sampling without replacement

import pandas as pd
import numpy as np

class immunoseq_graph:
    def __init__(self, immunoseq_data):
        self.immunoseq_data = immunoseq_data

    def subsample(self, immunoseq_data):
        t = []
        for row in df[sample].iterrows():
            t.extend([row[1].amino_acid] * row[1].templates)
        t = np.array(t)
        t = pd.DataFrame({
            'index': np.array(list(range(len(t)))),
            'amino_acid': t
        }).set_index('index')
        s_t = {}
        while len(s_t) < sampling_nr and len(t) > 0:
            s = t.sample(n=1, replace=False)
            if s.values[0][0] not in s_t:
                s_t[s.values[0][0]] = 1
                continue
            if s.values[0][0] in s_t:
                s_t[s.values[0][0]] += 1
            t = t.drop(s.index)










# Define input paths with sys args
sampling_nr = snakemake.input['']
in_file = snakemake.input['']
metadata = snakemake.input['']
public_sharing = snakemake.input['']
out_file = snakemake.output['']

# load AIRR seq data
df = pd.read_csv(in_file, sep='\t')

# get metadata from airr
nr_tcells = df.templates.sum()
nr_clones = len(df)

airr_meta_dict = {
    'sample': sample_col,
    'sample_nr': [int(sample.split('_')[1]) for sample in sample_col],
    'nr_tcells': nr_tcells_col,
    'nr_clones': nr_clones_col
}

# merge metadata and airr to one df
airr_meta = pd.DataFrame.from_dict(airr_meta_dict).sort_values(by='sample_nr',
                                                               ascending=True).drop('sample_nr',
                                                               axis=1).reset_index(drop=True)
airr_meta['clones_per_tcell'] = airr_meta.nr_clones / airr_meta.nr_tcells
airr_metadata = pd.merge(metadata, airr_meta, left_on='Sample', right_on='sample', ).drop('sample', axis=1)
airr_metadata.head()

# Sampling without replacement
df_sampled = {}

for sample in df.keys():
    print('sampling ', sample)
    t = []
    for row in df[sample].iterrows():
        t.extend([row[1].amino_acid] * row[1].templates)
    t = np.array(t)
    t = pd.DataFrame({
        'index': np.array(list(range(len(t)))),
        'amino_acid': t
    }).set_index('index')
    s_t = {}
    while len(s_t) < sampling_nr and len(t) > 0:
        s = t.sample(n=1, replace=False)
        if s.values[0][0] not in s_t:
            s_t[s.values[0][0]] = 1
            continue
        if s.values[0][0] in s_t:
            s_t[s.values[0][0]] += 1
        t = t.drop(s.index)

    s_t = pd.DataFrame.from_dict(s_t, orient='index').sort_values(by=0, ascending=False).reset_index()
    s_t.columns = ['amino_acid', 'templates']
    s_t['productive_frequency'] = s_t.templates / s_t.templates.sum()
    df_sampled[sample] = s_t
df_sampled_shared = {}
for dataset in df_sampled.keys():
    df_sampled_shared[dataset] = pd.merge(left=df_sampled[dataset], right=public_sharing,
                                          left_on='amino_acid', right_on='amino_acid', how='left').fillna(1)

from graph import RepertoireGraph

# create graphs
sample_col = []
ld_col = []
deg_dist_col = []
deg_freq_col = []
transi_col = []
authority_col = []
page_rank_col = []
eigen_centr_col = []
n_nodes_col = []
n_edges_col = []
max_deg_col = []
largest_comp_col = []
max_k_core_col = []
diameter_col = []
assortativity_col = []
average_degree_col = []
clust_coef_col = []
density_col = []
average_cluster_size_col = []

graph = {sample: RepertoireGraph(df_sampled_shared[sample]) for sample in df_sampled.keys()}
for s, g in graph.items():
    for ld in range(1, 2):
        # Make graph
        g.levenstein(ld=ld)
        g.threshold(ld=ld)
        g.make_graph(ld=ld)

        # append descriptivie statistics
        sample_col.append(s)
        ld_col.append(ld)
        deg_dist_col.append(g.degree_dist[ld])
        deg_freq_col.append(g.degree_freq[ld])
        transi_col.append(g.transitivity[ld])
        authority_col.append(g.authority[ld])
        page_rank_col.append(g.page_rank[ld])
        eigen_centr_col.append(g.eigen_centr[ld])
        n_nodes_col.append(g.n_nodes[ld])
        n_edges_col.append(g.n_edges[ld])
        max_deg_col.append(g.max_degree[ld])
        largest_comp_col.append(g.largest_comp[ld])
        max_k_core_col.append(g.max_k_core[ld])
        diameter_col.append(g.diameter[ld])
        assortativity_col.append(g.assortativity[ld])
        average_degree_col.append(g.average_degree[ld])
        clust_coef_col.append(g.clust_coef[ld])
        density_col.append(g.density[ld])
        average_cluster_size_col.append(g.average_cluster_size[ld])

# Merge into dict
desc_stats = {
    'sample': sample_col,
    'sample_nr': [int(nr.split('_')[1]) for nr in sample_col],
    'ld': ld_col,
    # 'deg_dist'       : deg_dist_col,
    # 'deg_freq'       : deg_freq_col,
    'transitivity': transi_col,
    'authority': authority_col,
    'page_rank': page_rank_col,
    'eigen_centr': eigen_centr_col,
    'n_nodes': n_nodes_col,
    'n_edges': n_edges_col,
    'max_degree': max_deg_col,
    'largest_comp': largest_comp_col,
    'max_k_core': max_k_core_col,
    'diameter': diameter_col,
    'assortativity': assortativity_col,
    'average_degree': average_degree_col,
    'clust_coef': clust_coef_col,
    'density': density_col,
    'average_cluster_size': average_cluster_size_col
}

# create dataframe with graph parameters and metadata
df_ds = pd.DataFrame.from_dict(desc_stats).sort_values(by=['ld', 'sample_nr']).drop('sample_nr', axis=1).reset_index(
    drop=True)
airr_metadata_graph = pd.merge(df_ds, airr_metadata, left_on='sample', right_on='Sample').drop('Sample', axis=1)
airr_metadata_graph.to_csv('/Users/eirikhoy/Dropbox/projects/airr_comet/data/graph_params_random.tsv', sep='\t',
                           index=False)