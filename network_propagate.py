import numpy as np 
from collections import defaultdict
import scipy.stats as stat
import pandas as pd
import networkx as nx 
from networkx.algorithms.link_analysis import pagerank
string_cutoff = 700
tmp_G = nx.Graph()
f1 = '9606.protein.info.v11.5.txt'
f2 = '9606.protein.links.v11.5.txt'
annotation = pd.read_csv(f1,header=0,sep='\t')
net = pd.read_csv(f2,header=0,sep=' ')
genemap=dict(zip(annotation.iloc[:,0],annotation.iloc[:,1]))
net['protein1']=net['protein1'].map(lambda x: genemap[x])
net['protein2']=net['protein2'].map(lambda x: genemap[x])
nodes1 = net.values[:,0]
nodes2 = net.values[:,1]
scores = net.values[:,2]
for n1, n2, score in zip(nodes1, nodes2, scores):
    if score >= string_cutoff:
        tmp_G.add_edge(n1, n2)
LCC_genes = max(nx.connected_components(tmp_G), key=len)
G = tmp_G.subgraph(LCC_genes)
network_nodes = G.nodes()
network_edges = G.edges()
#赋予权重
propagate_input = {}
drugtarget = []
for node in network_nodes:
    propagate_input[node] = 0
    if node in drugtarget:
        propagate_input[node] = 1
propagate_scores = pagerank(G, personalization=propagate_input, max_iter=100, tol=1e-06) ## NETWORK PROPAGATION
