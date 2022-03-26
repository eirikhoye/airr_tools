from igraph import *
import pandas as pd
import numpy as np
import stringdist
import imnet
import powerlaw
import igraph as ig
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import cm

def zero_to_one(n):
    if n == 0:
        return 1.0
    if n != 0:
        return n

class RepertoireGraph:

    def __init__(self, airr_df, aa_column = 'amino_acid', freq_column='productive_frequency'):
        self._junction_aa   = airr_df[aa_column].values
        self._freqs         = airr_df[freq_column].values
        self._sharing_level = list(map(zero_to_one, airr_df['sharing_level'])) # convert zero sharing level to 1
        self._dist_matrix   = {}
        self._bolean_matrix = {}
        self._graph         = {}
        self.degree_dist    = {}
        self.degree_freq    = {}
        self.transitivity   = {}
        self.authority      = {}
        self.page_rank      = {}
        self.eigen_centr    = {}
        self.n_nodes        = {}
        self.n_edges        = {}
        self.max_degree     = {}
        self.largest_comp   = {}
        self.max_k_core     = {}
        self.diameter       = {}
        self.assortativity  = {}
        self.average_degree = {}
        self.clust_coef     = {}
        self.density        = {}
        self.average_cluster_size = {}



    def levenstein(self, ld):
        """
        Create a levenstein distance matrix of input
        list of strings, for instance an immune receptor
        seqiemce repertoire (nt or aa, doesnt matter).
        """
        list_strings = self._junction_aa
        l = len(list_strings)
        dist_matrix = np.zeros((l, l))
        for idx, xval in enumerate(list_strings):
            for ydx, yval in enumerate(list_strings):
                d = stringdist.levenshtein(xval, yval)
                dist_matrix[idx, ydx] = d
        dist_matrix = pd.DataFrame(dist_matrix)
        dist_matrix.columns, dist_matrix.index = list_strings, list_strings
        self._dist_matrix[ld] = dist_matrix
        
    def threshold(self, ld):
        """
        Convert Levenstein Distance Matrix to binary 
        0 or 1 if distnance less than or equal to a 
        set levenstein threshiold (ld)
        0 means not connected clone 
        1 means connected clone
        """
        bolean_matrix = self._dist_matrix[ld]
        bolean_matrix[self._dist_matrix[ld] <= ld] = 1
        bolean_matrix[self._dist_matrix[ld] > ld] = 0
        self._bolean_matrix[ld] = bolean_matrix
        
    def make_graph(self, ld):
        """
        Make the iGraph from a given ld
        """
        self._graph[ld] = Graph.simplify(Graph.Adjacency((self._bolean_matrix[ld].values).tolist()))
        self._graph[ld].es['weight'] = self._bolean_matrix[ld].values[self._bolean_matrix[ld].values.nonzero()]
        self._graph[ld].vs['label']  = self._junction_aa
        self._graph[ld].vs['size_']   = self._freqs
        self._graph[ld].vs['sharing_level'] = self._sharing_level
        self._graph[ld].vs['degree'] = self._graph[ld].indegree()
        size_bol = np.array(self._graph[ld].vs['degree'])
        size_bol[size_bol >= 1] = 9e6
        size_bol[size_bol < 1] = 5e3
        self._graph[ld].vs['start_size'] = size_bol
        self._graph[ld].vs['size_'] = np.log(np.array(self._graph[ld].vs['size_']) * np.array(self._graph[ld].vs['start_size'])) + 2
        self._graph[ld].to_undirected()
        
        # Degree Distribution and Frequency
        self.degree_dist[ld] = sorted(self._graph[ld].vs['degree'], reverse=True)
        self.degree_freq[ld] = {i: self.degree_dist[ld].count(i) / len(self.degree_dist[ld]) for i in set(self.degree_dist[ld])}

        # Mean Local Parameters
        self.transitivity[ld] = self._graph[ld].transitivity_undirected()
        self.authority[ld]    = np.mean(self._graph[ld].authority_score())
        self.page_rank[ld]    = np.mean(self._graph[ld].personalized_pagerank(directed=False, damping=.85))
        self.eigen_centr[ld]  = np.mean(self._graph[ld].eigenvector_centrality(directed=False))

        # Global Parameters
        self.n_nodes[ld]        = len(self._graph[ld].vs)
        self.n_edges[ld]        = len(self._graph[ld].es)
        self.max_degree[ld]     = max(self.degree_dist[ld])
        self.largest_comp[ld]   = max(len(cl) for cl in self._graph[ld].clusters())
        self.max_k_core[ld]     = max(self._graph[ld].coreness())
        self.diameter[ld]       = self._graph[ld].diameter()
        self.assortativity[ld]  = self._graph[ld].assortativity_degree(directed=False)
        self.average_degree[ld] = np.mean(self._graph[ld].indegree())
        self.clust_coef[ld]     = np.mean(self._graph[ld].transitivity_local_undirected(mode='zero'))
        self.density[ld]        = np.mean(self._graph[ld].density(loops=False))
        self.average_cluster_size[ld] = np.mean([len(i) for i in self._graph[ld].clusters()])

    def plot_graph(self, ld):
        """
        Plot the graph, for a given ld
        """
        # Set up Visual style
        visual_style = {}    # Set bbox and margin
        visual_style['bbox'] = (800, 800)
        visual_style['margin'] = 25
        #visual_style["vertex_size"] = list(np.array(self._graph.vs['size']))    
        # Set vertex size
        visual_style["vertex_label_size"] = 0    # Set vertex lable size
        visual_style["edge_curved"] = False    # Don't curve the edges
        visual_style["edge_mode"] = '-'
        visual_style["edge_color"] = '#00000050'
        # Set vertex colors
        viridis = cm.get_cmap('viridis', int(np.max(np.unique(self._graph[ld].vs['sharing_level']))))
        newcolors = viridis(np.linspace(0, 1, int(np.max(np.unique(self._graph[ld].vs['sharing_level'])))))
        cmap = {int(i):list(newcolors[int(i)-1][0:3]) for i in np.unique(self._graph[ld].vs['sharing_level'])}
        #visual_style['vertex_color'] = [cmap[i] for i in self._graph[ld].vs['degree']]
        visual_style['vertex_color'] = [cmap[i] for i in self._graph[ld].vs['sharing_level']]
        # Set the Layout
        layout = self._graph[ld].layout_graphopt()
        visual_style['layout'] = layout
        y_coordinates = list(np.linspace(1, len(cmap.values()), len(cmap.values()))) # Added missing datapoint
        x_coordinates = list(np.repeat(3, len(cmap.values()))) # Added missing datapoint
        size_map      = list(np.repeat(300, len(cmap.values()))) # Added missing datapoint
        color_map = [color for color in list(cmap.values())[:len(x_coordinates)]]
        # Plot the graph!
        network_plot = {
            'network'  : ig.plot(self._graph[ld], **visual_style, inline=True),
            'colorbar' : plt.scatter(x_coordinates,y_coordinates, s = size_map, c = color_map)
        }
        return(network_plot)
        
        
    def powerlaw(self, ld, alt_model):
        results = powerlaw.Fit([i for i in self.degree_dist[[ld]] if i > 0])
        R, p = results.distribution_compare('power_law', alt_model)
        return(R, p)
        
        
        
# # plot barchart of top rearranged clones
#s = sample.head(10)
#chart = sns.barplot(data=s, x='junction_aa', y='freqs')
#for item in chart.get_xticklabels():
#    item.set_rotation(90)
#chart.set_title(str(Sample_ID))
#plt.ylim(0, 0.12)
#plt.plot()
#plt.show()

        