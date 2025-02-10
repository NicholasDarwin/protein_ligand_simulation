import networkx as nx
import matplotlib.pyplot as plt

def generate_interaction_map(interactions):
    G = nx.Graph()
    for interaction in interactions:
        G.add_edge(interaction[0], interaction[1], weight=interaction[2])
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True, node_size=2000)
    plt.show()
