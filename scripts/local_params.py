import pandas as pd
import numpy as np
import igraph as ig
import sys

filepath = sys.argv[1]
outpath = sys.argv[2]

################################################################################################

                           # Define function

################################################################################################

def local_params(g):
	local_param = {
	'transitivity' : np.mean(g.transitivity_local_undirected(mode='zero')),
	'authority' : np.mean(g.authority_score()),
	'pagerang' : np.mean(g.personalized_pagerank(directed=False)),
	'eigenvector_centrality' : np.mean(g.eigenvector_centrality(directed=False)),
	'closeness' : np.mean(g.closeness()),
	'betweenness' : np.mean(g.betweenness(directed=False))
	}

	local_param = pd.DataFrame(local_param, index=['parameter'])
	local_param = local_param.transpose()
	return local_param


#################################################################################################

                           # Load Graph

#################################################################################################

g = ig.Graph.Read_GraphML(filepath)

params = (local_params(g))
params.to_csv(sys.argv[2], sep='\t')
