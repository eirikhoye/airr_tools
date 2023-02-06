#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import os
import igraph as ig

# Load graphml files with igraph and store in dictionary
path_dir = 'data/imnet_output/graphs/'
graph_dict = {}
for file in os.listdir(path_dir):
    if not file.endswith('tsv.graphml'):
        continue
    sample_name = file.split('.')[0]
    g = ig.Graph.Read_GraphML(path_dir+file)
    graph_dict[sample_name] = g

# Create dataframe of global graph parameters, and store as global_param.tsv
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
global_param.to_csv(path_or_buf='data/global_param.tsv', sep='\t')

# Create dataframe of local graph parameters, and store as local_param.tsv
is_first=True
for dataset in graph_dict:

    local_param = {
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
local_param.to_csv(path_or_buf='data/local_param.tsv', sep='\t')

# Create dataframe for each clone with sharing level and degrees
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
degree_share.to_csv(path_or_buf='data/degree_share.tsv', sep='\t', index=False)
