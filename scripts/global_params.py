import pandas as pd
import numpy as np
import igraph as ig
import sys

filepath = sys.argv[1]
outpath = sys.argv[2]

################################################################################################

                           # Define function

################################################################################################

def global_params(g):
	global_param = {
	'nr_nodes' : len(g.vs),
	'nr_edges' : len(g.es),
	'largest_component' : max([len(i) for i in g.clusters()]),
	'maximal_k_core' : max(g.coreness()),
	'diameter' : g.diameter(),
	'assortativity' : g.assortativity_degree(directed=False),
	'average_cluster' : np.mean([len(i) for i in g.clusters()]),
	'average_degree' : np.mean(sorted(g.indegree(), reverse=True)),
	'clustering_coeff' : np.mean(g.transitivity_local_undirected(mode='zero')),
	'density' : np.mean(g.density(loops=False))
	}

	global_param = pd.DataFrame(global_param, index=['parameter'])
	global_param = global_param.transpose()
	return global_param


#################################################################################################

                           # Load Graph

#################################################################################################

g = ig.Graph.Read_GraphML(filepath)

params = (global_params(g))
params.to_csv(sys.argv[2], sep='\t')
