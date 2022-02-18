import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

G = nx.DiGraph()
edges = [ (332,343), (332,344), (333,332), (333,343), (333,349), (343,332), (343,333), (344,332), (344,333), (344,349), (349,332), (349,333), (349,344)]

G.add_edges_from(edges)
directed_edges = []
undirected_edges = []
for u, v in G.edges():
    if G.has_edge(v, u):
        undirected_edges.append((u, v))
    else:
        directed_edges.append((u, v))

pos = nx.spring_layout(G, iterations=50)
nx.draw_networkx_nodes(G, pos)
nx.draw_networkx_labels(G, pos)
nx.draw_networkx_edges(G, pos, edgelist=directed_edges)
nx.draw_networkx_edges(G, pos, edgelist=undirected_edges, arrows=False, edge_color='red')
plt.show()

