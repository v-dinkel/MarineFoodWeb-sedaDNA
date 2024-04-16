# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 14:28:10 2023

@author: vdinkel
"""

import numpy as np
import pandas as pd

def msum(A):
    return int(np.sum(np.sum(A)))

def binNorm(matrix):
    retM = matrix.where(matrix == 0, 1)
    return retM

def ESABO(x1, x2):
    # ESABO test
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
    # columns = taxa/features; rows = samples/observations
    
    ESABO_A = np.zeros((matrix.shape[1],matrix.shape[1])) # matrix[y][x]
    #ESABO_p = np.zeros((matrix.shape[0],matrix.shape[0])) # matrix[y][x]
    
    for i in range(0,matrix.shape[1]):
        for j in range(0,matrix.shape[1]):
            if i==j:
                ESABO_A[i][j] = 0
            else:
                arr1 = matrix[matrix.columns[i]].values # abundance vector of taxon i
                arr2 = matrix[matrix.columns[j]].values # abundance vector of taxon j
                if (0 in arr1 and 1 in arr1 and 0 in arr2 and 1 in arr2):
                    H_real, z_score = ESABO(arr1, arr2)
                    if np.isnan(z_score):
                        import pdb; pdb.set_trace()
                    ESABO_A[i][j] = z_score
                else:
                    #print ("not compatible with esabo")
                    ESABO_A[i][j] = 0
                
    ESABO_A = pd.DataFrame(ESABO_A)
    ESABO_A.columns = matrix.columns
    ESABO_A.index = matrix.columns
    
    ESABO_A = ESABO_A.where(abs(ESABO_A) > thresh, 0) # z score must be > 2
    #ESABO_A[ESABO_A != 0] = 1
    
    return ESABO_A#.astype(int)

def fillDiagonal(dataFrame):
    for index in dataFrame.columns:
        dataFrame.at[index,index] = 0
    return dataFrame

s_M = pd.read_csv("C:/Users/vdinkel/Desktop/Manuscript/input/F_KL-77.csv", delimiter=",", header=0, index_col=0)

#ESABO
s_M_bin = binNorm(s_M)
s_esaboM = corrESABOM(s_M_bin, 1.3)
s_esaboM = fillDiagonal(s_esaboM)
s_absEsabo = abs(s_esaboM.round().astype(int))
s_absEsabo = s_absEsabo.where(s_absEsabo == 0, 1)
print(msum(s_absEsabo))
pd.DataFrame.to_csv(s_absEsabo, "C:/Users/vdinkel/Desktop/Manuscript/output/KL-77_esabo.csv", sep=";")