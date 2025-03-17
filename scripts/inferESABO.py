# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 14:28:10 2023

@author: vdinkel
"""

import numpy as np
import pandas as pd

def msum(A):
    # returns the sum of a matrix
    return int(np.sum(np.sum(A)))

def binNorm(matrix):
    # normalizes the matrix to contain only 0 and 1 values (presence/absence)
    retM = matrix.where(matrix == 0, 1)
    return retM

def ESABO(x1, x2):
    # compute the ESABO score between two input vectors
    xij = np.logical_and(x1,x2).astype(int)
    p1 = np.sum(xij) / len(xij)
    
    p1a = np.sum(x1) / len(x1)
    p1b = np.sum(x2) / len(x2)
    
    std = np.sqrt(  (p1a * p1b * (1-(p1a*p1b))) / len(x1)  )
    mean = p1a * p1b

    z_score = (p1 - mean) / std
    
    H_real = mean
    return H_real, z_score

def corrESABOM(matrix, thresh):
    # compute the ESABO matrix given input abundance matrix and z-score threshold
    # columns = taxa/features; rows = samples/observations
    
    ESABO_A = np.zeros((matrix.shape[1],matrix.shape[1])) # matrix[y][x]
    
    for i in range(0,matrix.shape[1]):
        for j in range(0,matrix.shape[1]):
            if i==j:
                ESABO_A[i][j] = 0 # diagonal is 0
            else:
                arr1 = matrix[matrix.columns[i]].values # abundance vector of taxon i
                arr2 = matrix[matrix.columns[j]].values # abundance vector of taxon j
                if (0 in arr1 and 1 in arr1 and 0 in arr2 and 1 in arr2): # there must be presence/absence variation in the vectors
                    H_real, z_score = ESABO(arr1, arr2)
                    if np.isnan(z_score):
                        import pdb; pdb.set_trace()
                    ESABO_A[i][j] = z_score
                else:
                    #print ("not compatible with esabo")
                    ESABO_A[i][j] = 0
    
    # set the esabo indeces and columns to the column names of input matrix
    ESABO_A = pd.DataFrame(ESABO_A)
    ESABO_A.columns = matrix.columns
    ESABO_A.index = matrix.columns
    
    # apply threshold to ESABO adjacency matrix, set values below threshold to 0 
    ESABO_A = ESABO_A.where(abs(ESABO_A) > thresh, 0) # z score must be larger than given threshold
    
    return ESABO_A

def fillDiagonal(dataFrame):
    for index in dataFrame.columns:
        dataFrame.at[index,index] = 0
    return dataFrame

workdir ="/home/viktor/project_migration/MarineFoodWeb-sedaDNA/"
s_M = pd.read_csv(workdir+"input/F_KL-77.csv", delimiter=",", header=0, index_col=0)

#compute ESABO adjacency matrix and save into output folder
s_M_bin = binNorm(s_M)
s_esaboM = corrESABOM(s_M_bin, 1.3)
s_esaboM = fillDiagonal(s_esaboM)
s_absEsabo = abs(s_esaboM.round().astype(int))
s_absEsabo = s_absEsabo.where(s_absEsabo == 0, 1)
print(msum(s_absEsabo))
pd.DataFrame.to_csv(s_absEsabo, workdir+"output/KL-77_esabo.csv", sep=";")