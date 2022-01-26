# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from pyvis.network import Network
import networkx as nx
"""
Created on Wed Nov 17 11:55:11 2021

@author: nynkesimmes
"""


# NETWORK PYVIS - show graph

def networkxTest(df):
    net = Network(height='1000 px', width='1000 px', bgcolor='##34ebc3',
                  font_color='red', notebook=True,directed=True)
    net.barnes_hut()

    net.show_buttons(filter_=['physics'])
    net.from_nx(df['graphs'].values[0])
    net.show('graph_of_patient_test.html')

# networkxTest(df_small)

def show_graph(g, g_name):
    # net = Network(height = '1000 px', width = '1000 px', bgcolor='##34ebc3', font_color='red', notebook=True,directed=True)
    # net = Network(height ='1000 px', width ='1000 px', bgcolor='##34ebc3', font_color='red', notebook=True)
    # net.barnes_hut()
    # net.show_buttons(filter_=['physics'])
    # net.from_nx(g)
    # net.show(g_name)

    # nx.draw_networkx(g, with_labels = True)
    nx.draw(g, edge_vmin=10, edge_vmax=30, node_size=200, with_labels=True, font_size=7)
            # nodelist = ['Ataxia', 'Hypotonia', 'Delayed speech and language development', 'Motor delay']
    # plt.savefig(g_name)


def show_distribution(df):
    plt.hist(df['label'], rwidth=0.5)
    plt.title('Distribution of patients over affected genes')
    plt.xlabel('affected gene')
    plt.ylabel('frequency')
    plt.xticks(rotation=90)
    plt.show()


def show_similarity_matrix(measure_name, sim_matrix):
    fig, ax = plt.subplots(figsize=(20, 20))
    cax = ax.matshow(sim_matrix, interpolation='nearest')
    ax.grid(True)
    # plt.title('Similarity Matrix basd on', measure_name)
    plt.xticks(rotation=90)
    plt.yticks()
    fig.colorbar(cax, ticks=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, .75, .8, .85, .90, .95, 1])
    plt.show()    