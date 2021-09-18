#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 19:12:40 2020

@author: Oscar Rodriguez
"""
import numpy as np
import networkx as nx
from scipy import stats
import warnings
from itertools import combinations


warnings.filterwarnings('ignore');

def compute_spearman(data):
    '''
    Parameters
    ----------
    data : ndarray
        Array containig gene expression.
    name : string
        Name for the dataset. Used to identify output files.

    Returns
    -------
    rho : Sparse upper triangular csr_matrix
        Contains Spearman correlation coefficient for all samples.
        
    Computes the Spearman correlation for all samples.
    '''
    print("Computing Spearman correlation")
    rho, pval = stats.spearmanr(data, axis=1)
    rho = np.absolute(rho)
    rho = np.triu(rho, 1)
    print("Done!")
    return rho


def process_corr_matrix_using_adjacency_matrix(corr_matrix, threshold=0.2):
    '''

    Parameters
    ----------
    corr_matrix : ndarray
        correlation matrix
    threshold : float, optional
        Threshold to filter out connected genes. The default is 0.2.

    Returns
    -------
    corr_matrix : ndarray
        ndarray containing correlation coefficients above the threshold.
        
    This method computes an adjecency matrix from the correlation matrix
    and from this adjacency matrix it computes the largest connected component
    which is going to be the baseline network used downstream

    '''
    print("Processing correlation matrix:")
    adj_matrix = corr_matrix > threshold #computes de adjacency matrix
    G = nx.from_numpy_array(adj_matrix)
    #Now we extract all connected components from the adjacency matrix
    conn_comps = sorted(nx.connected_components(G), key=len, reverse=True)
    #We extract the largest connected component
    lcc = np.array(list(conn_comps.pop(0)))
    print("Largest connected component: "+str(len(lcc)))
    # We get rid of the entries on the adjacency matrix that correspond 
    # to the smaller connected components
    for index_set in conn_comps:
        if len(index_set)>1:
            idxs = combinations(index_set,2)
            for i_j in idxs:
                adj_matrix[i_j] = False
    # Now adjacency_matrix is a mask, we use it to filter out all genes with
    # correlation below the threshold using entry-wise multiplication
    corr_matrix = corr_matrix*adj_matrix
    print("Done!")
    return corr_matrix