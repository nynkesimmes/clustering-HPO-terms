# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
"""
Created on Mon Apr 14 09:01:18 2014

"""


def clusterPlot(X, clusterid, centroids=None, y=None):
    '''
    CLUSTERPLOT Plots a clustering of a data set as well as the true class
    labels. If data is more than 2-dimensional it should be first projected
    onto the first two principal components. Data objects are plotted as a dot
    with a circle around. The color of the dot indicates the true class,
    and the cicle indicates the cluster index. Optionally, the centroids are
    plotted as filled-star markers, and ellipsoids corresponding to covariance
    matrices (e.g. for gaussian mixture models).

    Usage:
    clusterplot(X, clusterid)
    clusterplot(X, clusterid, centroids=c_matrix, y=y_matrix)
    clusterplot(X, clusterid, centroids=c_matrix, y=y_matrix, covars=c_tensor)

    Input:
    X           N-by-M data matrix (N data objects with M attributes)
    clusterid   N-by-1 vector of cluster indices
    centroids   K-by-M matrix of cluster centroids (optional)
    y           N-by-1 vector of true class labels (optional)
    '''

    X = np.asarray(X)
    cls = np.asarray(clusterid)
    if y is None:
        y = np.zeros((X.shape[0], 1))
    else:
        y = np.asarray(y)
    if centroids is not None:
        centroids = np.asarray(centroids)
    K = np.size(np.unique(cls))
    C = np.size(np.unique(y))
    ncolors = np.max([C, K])

    # CREATE labels_true_int, dus true class labels in vorm van integers
    current_int = -1
    y_int = [0]*len(y)
    for i, label in enumerate(y):
        if label != y[i-1]:
            current_int += 1
        y_int[i] = current_int
    # print(y_int)

    # VERSION WITH COLORS
    colors = [0]*ncolors
    for color in range(ncolors):
        colors[color] = plt.cm.jet.__call__((color*255)//(ncolors-1))[:3]
    for i, cs in enumerate(np.unique(y)):
        plt.plot(X[(y == cs).ravel(), 0], X[(y == cs).ravel(), 1], 'o',
                 markeredgecolor='k', markerfacecolor=colors[i], markersize=6,
                 zorder=2)
    for i, cr in enumerate(np.unique(cls)):
        plt.plot(X[(cls == cr).ravel(), 0], X[(cls == cr).ravel(), 1], 'o',
                 markersize=11, markeredgecolor=colors[i],
                 markerfacecolor='None', markeredgewidth=2, zorder=1)

    # VERSION WITH MARKERS
    # colors = [0]*ncolors
    # markers = ['o', '.', ',', 'x', '+', 'v', '^', '<', '>', 's', 'd', '*', 'D', 'P', '8', '_']

    # for color in range(ncolors):
    #     colors[color] = plt.cm.jet.__call__((color*255)//(ncolors-1))[:3]

    # y_markers = [markers[i] for i in list(np.array(y_int)%16)]
    # cls_colors = [colors[i] for i in list(np.array(clusterid)%16)]
    # # print(y_markers)
    # # print(cls_colors)

    # for i in range(X.shape[0]):
    #     plt.scatter(X[i,0], X[i,1], color=cls_colors[i], marker=y_markers[i],
    #                 edgecolors='k', alpha=1, s=mpl.rcParams['lines.markersize'] ** 2.2)

    plt.xlabel('1-eigenvector')
    plt.ylabel('2-eigenvector')

    # create legend
    legend_items = (np.unique(y).tolist() + np.unique(cls).tolist() +
                    np.unique(cls).tolist())
    for i in range(len(legend_items)):
        if i < C:
            legend_items[i] = 'Class: {0}'.format(legend_items[i])
        elif i < C + K:
            legend_items[i] = 'Cluster: {0}'.format(legend_items[i])
        else:
            legend_items[i] = 'Centroid: {0}'.format(legend_items[i])
    plt.legend(legend_items, numpoints=1, markerscale=1, prop={'size': 9},
               ncol=2, bbox_to_anchor=(1.05, 1), loc='upper left')
