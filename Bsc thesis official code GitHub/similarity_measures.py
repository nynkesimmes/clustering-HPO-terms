# -*- coding: utf-8 -*-
import numpy as np
import networkx as nx
import gmatch4py as gm
"""
Created on Wed Nov 17 11:35:54 2021

@author: nynkesimmes
"""


def rescale_matrix(matrix):
    """
    Rescale/normalize matrix to values between 0 and 1.

    Parameters
    ----------
    matrix : Array of float64
        Matrix with float values.

    Returns
    -------
    matrix_rescaled : Array of float64
        Matrix with float values between 0 and 1.

    """
    max_val = np.max(matrix)
    min_val = np.min(matrix)
    matrix_rescaled = np.array([(x-min_val)/(max_val-min_val) for x in matrix])

    return matrix_rescaled


def create_sem_sim_matrix(df, scorer):
    """
    Create similarity matrix based on semantic similarity measure.

    This function calculates the semantic similarity between every data point
    with the scorer that is given as parameter e.g. Resnik or Lin. Score is
    saved in a matrix.

    Parameters
    ----------
    df : DataFrame
        Dataframe of 553 patients with HPO terms and corresponding HPO graphs.
    scorer : Scorer
        Instantiance of class scorer assigned with specific measure.

    Returns
    -------
    sim_matrix : Array of float64
        Similarity matrix resulting from Resnik or Lin.

    """
    sim_matrix = np.zeros([len(df), len(df)])
    for i, terms1 in enumerate(df['hpo_all']):
        for j, terms2 in enumerate(df['hpo_all']):
            sim_matrix[i][j] = scorer.score_term_sets_basic(terms1, terms2)

    if scorer.scoring_method == 'Resnik':  # rescale to values between 0 and 1
        sim_matrix = rescale_matrix(sim_matrix)

    return sim_matrix


def create_graph_sim_matrix_ged_gm(df):
    """
    Create similarity matrix based on graph edit distance.

    This function uses method from library GMatch4py to calculated the GED
    between a whole list of graphs.

    Parameters
    ----------
    df : DataFrame
        Dataframe of 553 patients with HPO terms and corresponding HPO graphs.

    Returns
    -------
    graph_sim_matrix_GED : Array of float64
        Similarity matrix resulting from GED measure.

    """
    print("Start creating graph similarity matrix based on GED")
    ged = gm.GraphEditDistance(1, 1, 0, 0)  # node edit costs are 1, edge edit costs are 0
    graph_sim_matrix_GED = ged.similarity(ged.compare(list(df['graphs']), None))  # pass all patient graphs at once

    return graph_sim_matrix_GED


def getMCS(g1, g2):
    """
    Get Maximum Common Subgraph of two given graphs.

    This function checks for every edge of one graph whether it also occurs in
    the other graph. Its induced subgraph togheter with its corresponding size
    is returned.

    Parameters
    ----------
    g1 : graph
        First graph.
    g2 : graph
        Second graph.

    Returns
    -------
    (graph, int)
        The induced subgraph of g1 and g2 and the largest component.

    """
    matching_graph = nx.Graph()

    for n1, n2 in g2.edges():
        if g1.has_edge(n1, n2):
            matching_graph.add_edge(n1, n2)

    largest_component = max(nx.connected_components(matching_graph), key=len)
    return nx.induced_subgraph(matching_graph, largest_component)


def create_graph_sim_matrix_mcs(df):
    """
    Create similarity matrix based on Maximum Common Subgraph.

    This function gets a data set also containing graphs which are all compared
    with each other by use of helper function getMCS. A similarity matrix is
    created based on the relative number of nodes in the MCS of every two
    graphs.

    Parameters
    ----------
    df : DataFrame
        Dataframe of 553 patients with HPO terms and corresponding HPO graphs.

    Returns
    -------
    graph_sim_matrix_MCS_nnodes_relative_avg : Array of float65
        Similarity matrix resulting from MCS measure.

    """
    print("Start creating graph similarity matrix based on MCS")
    graph_sim_matrix_MCS_nnodes_relative_avg = np.zeros([len(df), len(df)])

    nr_nodes_list = []
    for i, graph1 in enumerate(df['graphs']):
        for j, graph2 in enumerate(df['graphs']):
            min_nnodes12 = min(graph1.number_of_nodes(), graph2.number_of_nodes())
            max_nnodes12 = max(graph1.number_of_nodes(), graph2.number_of_nodes())
            avg_nnodes12 = (min_nnodes12 + max_nnodes12) / 2
            mcs12 = getMCS(graph1, graph2)
            mcs12_nnodes = mcs12.number_of_nodes()
            nr_nodes_list.append(mcs12_nnodes)

            graph_sim_matrix_MCS_nnodes_relative_avg[i][j] = mcs12_nnodes/avg_nnodes12

    return graph_sim_matrix_MCS_nnodes_relative_avg
