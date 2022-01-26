# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import SpectralClustering
from clusterPlot import clusterPlot
from sklearn.metrics import adjusted_rand_score
from sklearn.cluster import KMeans
import numpy as np

"""
Created on Wed Nov 17 11:52:04 2021

@author: nynkesimmes
"""

def calculate_ARI(measure_name, labels_true, labels_pred):
    """
    Calculate ARI of the true and predicted clustering result.

    This function calculated the Adjusted Rand Index on the basis of the true
    clustering result (labels_true) and the predicted clustering result
    (labels_pred).

    Parameters
    ----------
    measure_name : str
        String that indicates the concerned similarity measure.
    labels_true : Array of object
        Array that contains all the true labels i.e. the actual syndromes.
    labels_pred : Array of int64
        Array that contains the predicted cluster numbers.

    Returns
    -------
    ARI : float64
        The Adjusted Rand Index of the true and predicted labels.

    """
    ARI = adjusted_rand_score(labels_true, labels_pred)
    # print("Rand index of", measure_name, RI)
    print("ARI of ", measure_name, ARI)
    return ARI


def cluster_spectral(measure_name, sim_matrix, labels_true, nr_labels):
    """
    Apply spectral clustering to given similarity matrix.

    Parameters
    ----------
    measure_name : str
        String that indicates the concerned similarity measure..
    sim_matrix : Array of float64
        Similarity matrix that is given as parameter.
    labels_true : Array of object
        Array that contains all the true labels i.e. the actual syndromes.
    nr_labels : int
        Integer that indicates the number of classes / clusters the original
        data set consists of.

    Returns
    -------
    ARI : float64
        The Adjusted Rand Index of the true and predicted labels.

    """
    # Apply pca on 2 principal components
    pca = PCA(2)
    sim_matrix_pca2 = pca.fit_transform(sim_matrix)
    spectral = SpectralClustering(n_clusters=nr_labels,
                                  assign_labels='discretize',
                                  random_state=0,
                                  affinity='precomputed')
    spectral = spectral.fit(sim_matrix)  # apply spectral clustering to the whole matrix
    labels_pred_spectral = spectral.fit_predict(sim_matrix)  # get the predicted labels i.e. cluster ids

    # Calculate ARI
    ARI = calculate_ARI(measure_name, labels_true, labels_pred_spectral)

    # Create cluster plot
    clusterPlot(sim_matrix_pca2, labels_pred_spectral, y=labels_true)
    plt.title(measure_name + " with ARI " + str(round(ARI, 2)))
    plt.show()

    return ARI


def cluster_kmeans(measure_name, sim_matrix, labels_true, nr_labels):
    """
    Apply kmeans clustering to given similarity matrix.

    Parameters
    ----------
    measure_name : str
        String that indicates the concerned similarity measure..
    sim_matrix : Array of float64
        Similarity matrix that is given as parameter.
    labels_true : Array of object
        Array that contains all the true labels i.e. the actual syndromes.
    nr_labels : int
        Integer that indicates the number of classes / clusters the original
        data set consists of.

    Returns
    -------
    ARI : float64
        The Adjusted Rand Index of the true and predicted labels.

    """
    kmeans = KMeans(n_clusters=nr_labels).fit(sim_matrix)
    labels_pred_kmeans = kmeans.predict(sim_matrix)
    centroids = kmeans.cluster_centers_

    # Calculate ARI
    ARI = calculate_ARI(measure_name, labels_true, labels_pred_kmeans)

    # Create cluster plot
    clusterPlot(sim_matrix, labels_pred_kmeans, centroids, labels_true)
    plt.title(measure_name + " with ARI " + str(round(ARI, 2)))
    plt.show()

    return ARI


# PRINT TRUE EN PREDICTED LABELS NEXT TO EACH OTHER
def compare_labels(measure_name, labels_true, labels_pred):
    """
    Print labels of original and predicted clustering result.

    Parameters
    ----------
    measure_name : str
        String that indicates the concerned similarity measure.
    labels_true : Array of object
        Array that contains all the true labels i.e. the actual syndromes.
    labels_pred : Array of int64
        Array that contains the predicted cluster numbers.

    Returns
    -------
    None.

    """
    current_int = 0
    labels_true_int = np.zeros((len(labels_true),), dtype=int)
    for i, label in enumerate(labels_true):
        if label != labels_true[i-1]:
            current_int += 1
        labels_true_int[i] = current_int

    for i in range(len(labels_true)):
        print(labels_true[i], labels_true_int[i], labels_pred[i])
