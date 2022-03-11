from click import option
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from pyvis.network import Network

G = nx.DiGraph()
f = open("scripts/plot_data", "r")
data = f.read().split('\n')
edges = []
for x in data:
    xs = x.split(" ")
    if len(xs) == 2:
        edges.append((int(xs[0]), int(xs[1])))

G.add_edges_from(edges)

directed_edges = []
undirected_edges = []
for u, v in G.edges():
    if G.has_edge(v, u):
        undirected_edges.append((u, v))
    else:
        directed_edges.append((u, v))

G.clear()
G.add_edges_from(directed_edges)
G.add_edges_from(undirected_edges, color='red', arrows='false')

net = Network('720px', '1280px', notebook=True, directed=True)

net.from_nx(G)
net.show_buttons(filter_=['physics'])
net.show("example.html")

# pos = nx.spring_layout(G, iterations=50)
# nx.draw_networkx_nodes(G, pos)
# nx.draw_networkx_labels(G, pos)
# nx.draw_networkx_edges(G, pos, edgelist=directed_edges)
# nx.draw_networkx_edges(G, pos, edgelist=undirected_edges, arrows=False, edge_color='red')
# plt.show()

