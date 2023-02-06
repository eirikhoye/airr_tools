import findspark
findspark.init()
import pyspark
import os
import sys
import imnet
import numpy as np
import pandas as pd
import networkx as nx

"""
Usage: python scripts/run_imnet.py {path/to/input.tsv} {name_outfile}
"""


os.environ['SPARK_HOME'] = '/opt/spark/'
sc = pyspark.SparkContext('local[*]')

def run_imnet(file_path, out_path, min_ld=1, max_ld=1):
    df = pd.read_csv(file_path, sep='\t')
    strings = list(df.amino_acid)

    g = imnet.process_strings.generate_graph(strings, sc=sc, min_ld=1, max_ld=1)
    g.remove_edges_from(nx.selfloop_edges(g))   
    param_glob = {}

    nx.write_graphml(g, out_path+'.graphml')
    param_glob['Nodes(V)'] = g.number_of_nodes()
    param_glob['Edges(E)'] = g.number_of_edges()
    param_glob['LargestComponent'] = len(max(nx.connected_components(g), key=len))
    param_glob['max_k_core'] = len(nx.k_core(g, k=None, core_number=None).edges())
    param_glob['clique'] = nx.graph_clique_number(g)
    print(param_glob)

    degrees = imnet.process_strings.generate_degrees(strings, sc=sc, min_ld=min_ld, max_ld=max_ld)
    np.savetxt(out_path+'_degree.txt', degrees, fmt='%1i')



    


run_imnet(sys.argv[1], sys.argv[2], min_ld=1, max_ld=1)


sc.stop()


