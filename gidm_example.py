#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 23:27:19 2021

@author: Oscar Rodriguez
"""
from scipy.io import loadmat
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import save_npz, load_npz
import multiprocessing as mp
import itertools
import time
import os
from gidmutils import *
import gc

fpath = '../TestData/Cancer/outputmat/'
os.makedirs(os.path.dirname(fpath), exist_ok=True)

def fill_sparse(i):
    try:
        sample=data[i,:]
        gene_idx_1=[]
        gene_idx_2=[]
        vals=[]
        nz=np.nonzero(sample)
        nz=nz[0]
        for k in range(len(nz)-1):
            t=len(nz)-1 - k
            temp=[[nz[k]]*t]
            temp2=nz[k+1:len(nz)]
            gene_idx_2.append(temp2)
            temp=list(itertools.chain(*temp))
            gene_idx_1.append(temp)
            vals.append(0.5*rho[temp, temp2]*(sample[temp]+sample[temp2]))  
        gene_idx_1=list(itertools.chain(*gene_idx_1))
        gene_idx_2=list(itertools.chain(*gene_idx_2))
        vals=list(itertools.chain(*vals))
        u=csr_matrix((vals, (gene_idx_1, gene_idx_2)), shape=(data.shape[1], data.shape[1]), dtype=np.float32)
        filename=fpath+str(i)+'.npz'
        save_npz(filename, u)
        return 1
    except Exception as e:
        print(e)
        
def get_distance(i, j):
    try:
        ci=load_npz(fpath+str(i)+'.npz')
        cj=load_npz(fpath+str(j)+'.npz')
        diff = ci-cj
        diff_v = np.abs(diff.data)
        d=np.sum(diff_v)
        return [i,j,d]
    except Exception as e:
        print(e)


scData_m = loadmat('../TestData/Cancer/breast_cancer_5000.mat')
data = scData_m['A']

rho = compute_spearman(data.T)
rho = process_corr_matrix_using_adjacency_matrix(rho, 0.37)

try:
    #The next line is to set the number of CPUs
    #Change it accordingly if you are using a SLURM managed cluster
    N_CPUs =4#int(os.getenv('SLURM_CPUS_ON_NODE'))
    
    pool = mp.Pool(N_CPUs)
    result_objects = []
    print("Computing sparse matrices...")
    tic = time.perf_counter()
    for i in range(data.shape[0]):
        result_objects.append(pool.apply_async(fill_sparse, args=(i,)))
    pool.close()
    pool.join()
    toc = time.perf_counter()    
    results = [r.get() for r in result_objects]
    results=[]
    result_objects=[]
    print(f"Done: {toc - tic:0.4f} seconds")
    del rho
    gc.collect()
    pool = mp.Pool(N_CPUs)
    tic = time.perf_counter()
    print("Computing distances...")
    for i in range(data.shape[0]-1):
        for j in range(i+1, data.shape[0]):
            result_objects.append(pool.apply_async(get_distance, args=(i, j)))
    pool.close()
    pool.join()
    results = [r.get() for r in result_objects]
    D = np.zeros((data.shape[0], data.shape[0]))
    for r in results:
        D[r[0], r[1]] = r[2]
    toc = time.perf_counter()
    print(f"Distances done: {toc - tic:0.4f} seconds")
except Exception as e:
    print(e)

print("Saving distance matrix as cancerMP_21_signed")
np.save('../TestData/Cancer/cancerMP_21_signed.npy', D)
print("Cancer completed!")