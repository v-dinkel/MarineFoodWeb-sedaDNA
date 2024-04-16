# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 13:38:17 2022

@author: vdinkel
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx
import networkx.algorithms.community as nx_comm
import scipy.stats
from collections import Counter
import statsmodels.formula.api as smf


def msum(A):
    return int(np.sum(np.sum(A)))

def AND(A, B):
    return B.where(A == 1, 0)

def OR(A, B):
    C = abs(A) + abs(B)
    return C.replace(2, 1)

def fillDiagonal(dataFrame):
    for index in dataFrame.columns:
        dataFrame.at[index,index] = 0
    return dataFrame

def getOwnTPs(trueA, predA, threshold):
    tp = 0;fp = 0; tn = 0; fn = 0
    
    i = 0
    for predArr in predA.values:
        j = 0
        for predVal in predArr:
            predVal = abs(predVal)
            if predVal < threshold:
                predVal = 0
            if predVal > 0 and trueA.values[i][j] > 0:
                tp += 1
            if predVal > 0 and trueA.values[i][j] == 0:
                fp += 1
            if predVal == 0 and trueA.values[i][j] == 0:
                tn += 1
            if predVal == 0 and trueA.values[i][j] > 0:
                fn += 1
            j+=1
        i += 1
    if tp+fn != 0:
        tpr = tp/(tp+fn)
    else:
        tpr = 0
    if fp+tn != 0:
        fpr = fp/(fp+tn)
    else:
        fpr = 0
    if tp+fp+tn+fp != 0:
        acc = (tp+tn) / (tp+fp+tn+fp)
    else:
        acc = 0
    return {"tp": tp, "fp":fp, "tn": tn, "fn": fn, "tpr": tpr, "fpr": fpr, "acc": acc}

def compareNetworks(baseA, testB):
    intersect_all = AND(abs(baseA), abs(testB))
    try:
        JI_all = (msum(intersect_all)) / min(msum(abs(baseA)), msum(abs(testB))) #(msum(abs(baseA)) + msum(abs(testB)) - msum(abs(intersect_all)))
    except:
        JI_all = 0
    return JI_all

def getNetworkStats(M):
    G = nx.from_numpy_matrix(np.matrix(M))
    print ("Stats")
    
def modularity(M):
    if type(M) == type(pd.DataFrame()):
        G = nx.from_numpy_matrix(np.matrix(M))
    else:
        G = M
    comms = nx_comm.greedy_modularity_communities(G, resolution=1.0)
    return float(round(nx_comm.modularity(G, comms),2)), comms

def cc(M):
    if type(M) == type(pd.DataFrame()):
        G = nx.from_numpy_matrix(np.matrix(M))
    else:
        G = M
    return float(round(nx.average_clustering(G), 2))

def plotG(M):
    plt.clf()
    G = nx.from_numpy_matrix(np.matrix(M))
    nx.draw(G)
    return

def calc_z(x, shuffled):
    if np.std(shuffled) != 0:
        return round((x - np.mean(shuffled)) / np.std(shuffled),2)
    else:
        return 0.00
    
def getNodeModuleID(node, modules):
    for i in range(0, len(modules['comms'])):
        if node in modules['comms'][i]:
            return i
            
def exportToGephi(G, df_families, tpos, tneg, modules, fileName, trophicLevels):
    f = open(fileName, "w+")
    f.write("graph\n[\n")
    # NODES
    
    degrees = dict(G.degree)
    nodes = list(set([k if degrees[k] > 0 else None for k in degrees.keys()]))
    for node in nodes:
        if node != None:
            moduleID = getNodeModuleID(node, modules)
            f.write("node\n[\n")
            f.write('id "'+node+'"\n')
            f.write('label "'+node+'"\n')
            f.write('module "'+str(moduleID)+'"\n')
            #f.write('abund '+str(int(np.sqrt(np.sum(abund[node]))))+'\n')
            #f.write('zeros '+str(int(analysis_dict['nzeros'][node]))+'\n')
            
            try:
                taxGroup = (df_families[df_families[0] == node][1]).values[0].replace("heterotrophic protitsts", "heterotrophic protists").replace("crustaceous zooplnakton", "crustaceous zooplankton")
                f.write('taxonomic_group "'+taxGroup+'"\n')
                f.write('trophic_level "'+str(trophicLevels[trophicLevels[0] == taxGroup][1].values[0])+'"\n')
            except:
                #import pdb; pdb.set_trace()
                f.write('taxonomic_group "NaN"\n')
                f.write('trophic_level "NaN"\n')
            
            # SCORRS
            scorr = 0
            sval = 0
            if node in tpos:
                f.write('scorr 1\n')
            elif node in tneg:
                f.write('scorr -1\n')
            else:
                f.write('scorr 0\n')
            
            f.write("]\n")
    
    # EDGES
    for edge in G.edges:
        #scndIndex = np.where(pp_cgr_A.columns == edge[1])[0][0]
        #if pp_cgr_A[edge[0]][scndIndex] == -1:
        #    import pdb; pdb.set_trace()

        #import pdb; pdb.set_trace()
        #weight = pp_cgr_A_edgelist[(pp_cgr_A_edgelist['from'] == edge[0]) & (pp_cgr_A_edgelist['to'] == edge[1])]['weight'].values[0]
        #weight = G.get_edge_data(edge[0],edge[1])['weight']
        
        f.write("edge\n[\n")
        f.write('source "'+edge[0]+'"\n')
        f.write('target "'+edge[1]+'"\n')
        

        f.write('weight '+str(G.get_edge_data(edge[0],edge[1])['weight'])+'\n')

            
        #if weight>=0:
        #    f.write('assoc 1\n')
        #else:
        #    f.write('assoc -1\n')
        
        #f.write('intersect_all '+str(int(analysis_dict['intersect_all'][edge[0]][edge[1]]))+'\n')
        #f.write('intersect_pos '+str(int(analysis_dict['intersect_pos'][edge[0]][edge[1]]))+'\n')
        #f.write('intersect_neg '+str(int(analysis_dict['intersect_neg'][edge[0]][edge[1]]))+'\n')
        #import pdb; pdb.set_trace()

        f.write("]\n")
    
    f.write("]")
    f.close()

def modcorrs(M, poscorrs, negcorrs, thresh):
    G = nx.from_pandas_adjacency(M)
    
    mod, test_comms = modularity(G)
    com_lengths = [len(k) for k in list(test_comms)]
    sorted_comm_ind = sorted( com_lengths, reverse=True)
    checked_comsizes = []
    accoccs = 0
    badoccs = 0
    distcomms = 0
    mixedcomms = 0
    for k in sorted_comm_ind:
        if k > 1 and k not in checked_comsizes:
            indeces = [i for i, x in enumerate(com_lengths) if x == k]
            for ind in indeces:
                thiscom = list(list(test_comms)[ind])
                posocc = len(set(thiscom) & set(poscorrs))
                negocc = len(set(thiscom) & set(negcorrs))
                
                if (negocc + posocc > 0):
                    fract = 0.0
                    if max(posocc, negocc) > 0:
                        fract = min(posocc, negocc) / max(posocc, negocc)
                    if fract == 0.0:
                        #import pdb; pdb.set_trace()
                        accoccs += posocc + negocc
                        if (posocc + negocc) > 1 and len(thiscom) >= 3 :
                            distcomms += 1
                    else:
                        badoccs += posocc + negocc
                        mixedcomms += 1
                    #print ("n: ", len(thiscom), " | posocc: ", posocc, " | negocc: ", negocc, " | fract: ", round(fract,2))
            checked_comsizes.append(k)
    goodresult = False
    if (distcomms > mixedcomms) and (accoccs > badoccs) and (distcomms - mixedcomms)>=3 and ((badoccs == 0 ) or (accoccs/badoccs)>=thresh):
        goodresult = True
    return {"nodes": len(M.columns),"edges": int(msum(M)/2),"modularity": mod, "cc":cc(M), "comms": test_comms, "distcomms": distcomms, "mixedcomms": mixedcomms, "accoccs": accoccs, "badoccs": badoccs, "goodresult": goodresult }


def getModuleMTL(module):
    #getModuleMTL(mods2["comms"][0])
    this_module = list(module)
    mtl = int(np.mean([trDict2[k] for k in this_module])*100)/100
    
    relTrSums = []
    for trLvl in trDict.keys():
        trSet = list(set(this_module).intersection(trDict[trLvl]))
        relTrSums.append(int(abunds[trSet].sum(axis=1).sum(axis=0) *10000)/10000 )
    print (mtl, relTrSums)
    return mtl, relTrSums

def getModuleTemp(module):
    thisMod = list(module)
    positives = list(set(thisMod).intersection(tpos))
    negatives = list(set(thisMod).intersection(tneg))
    print ( {"pos": len(positives), "neg": len(negatives)})
    return {"pos": len(positives), "neg": len(negatives)}

def getGlobiMatrix(df_globi, columns, taxa, rank):
    df_sub_globi = df_globi.loc[df_globi['s_'+rank].isin(taxa) & df_globi['t_'+rank].isin(taxa)]
    M_globi = pd.DataFrame(columns=sorted(columns), index = sorted(columns)).fillna(0)
    df_sub_globi.index = [k for k in range(0, len(df_sub_globi))]
    for idx in df_sub_globi.index:
        M_globi.at[df_sub_globi.loc[idx]['s_'+rank], df_sub_globi.loc[idx]['t_'+rank]] = 1
        M_globi.at[df_sub_globi.loc[idx]['t_'+rank], df_sub_globi.loc[idx]['s_'+rank]] = 1
    return M_globi

def shuffleNetwork(M, nshuffles, temps, globi):
    modularities = []
    ccs = []
    
    icover_jis = []
    icover_TPs = []
    icover_FPs = []
    icover_FNs = []
    
    distinct_comms = []
    mixed_comms = []
    distinct_taxa = []
    mixed_taxa = []
    goodresults = []
    
    for i in range(0, nshuffles):
        G = nx.from_pandas_adjacency(M)
        G_shuffled = nx.double_edge_swap(G, nswap=len(G.edges)*10, max_tries=len(G.edges)*1000)
        M_shuffled = nx.adjacency_matrix(G).todense()
        M_shuffled = pd.DataFrame(M_shuffled, columns = M.columns, index = M.index)
        
        ccs.append(cc(G_shuffled))
        
        iCover_ji = compareNetworks(M_shuffled, globi)
        icover_jis.append(iCover_ji)
        
        iCover_Own = getOwnTPs(globi, M_shuffled, 0.1)
        icover_TPs.append(iCover_Own["tp"])
        icover_FPs.append(iCover_Own["fp"])
        icover_FNs.append(iCover_Own["fn"])
        
    return {"modularities": modularities, "css": ccs, "icover_jis": icover_jis, "icover_TPs": icover_TPs, "icover_FPs": icover_FPs, "icover_FNs": icover_FNs, "distcomms":distinct_comms, "mixedcomms" : mixed_comms, "distoccs" : distinct_taxa, "mixedoccs" : mixed_taxa, "goodresults" : goodresults}

def getMultiplexThresh(multA, thresh):
    multA = multA.where(multA >= thresh, 0)
    multA = multA.where(multA == 0, 1)
    return multA

def AMI_compartments(G, ascDicts, vabunds):
    # T(i+) is the sum of all flows leaving species i, and T(+j) is the sum of all flows entering species j.
    A = 0
    T = len(list(G.edges))
    ascEdges = ascDicts[0]
    nAscEdges = ascDicts[1]
    ascEdges_list = ascDicts[2]
    
    degs = G.degree()
    vabunds = vabunds * 1
    connectedNodes = []
    for edge in list(G.edges):
        connectedNodes.append(edge[0])
        connectedNodes.append(edge[1])
        
        if vabunds[edge[0]].values[0] == 0.0 or vabunds[edge[1]].values[0] == 0.0:
            G[edge[0]][edge[1]]["weight"] = 0.0
        else:
            e0_zeroNeighbour = sum([vabunds[e[0]].values[0] == 0.0 or vabunds[e[1]].values[0] == 0.0 for e in G.edges(edge[0])])
            e1_zeroNeighbour = sum([vabunds[e[0]].values[0] == 0.0 or vabunds[e[1]].values[0] == 0.0 for e in G.edges(edge[1])])
            e0 = vabunds[[edge[0]]].values[0][0] / (degs[edge[0]] - e0_zeroNeighbour)
            e1 = vabunds[[edge[1]]].values[0][0] / (degs[edge[1]] - e1_zeroNeighbour)
            G[edge[0]][edge[1]]["weight"] = e0+e1
            
    connectedNodes = list(set(connectedNodes))
    T = np.sum(np.sum(vabunds[connectedNodes]))
    
    for key in ascEdges.keys():
        ascEdges[key] = np.sum([G[edge[0]][edge[1]]["weight"] for edge in ascEdges_list[key]])
    #import pdb; pdb.set_trace()
    
    try:
        T12 = ascEdges["1_to_2"]/T * np.log2( (ascEdges["1_to_2"]*T)/((ascEdges["1_to_2"]+ascEdges["1_to_3"]+ascEdges["1_to_4"])*(ascEdges["1_to_2"])))
        H12 = ascEdges["1_to_2"]/T * np.log2(ascEdges["1_to_2"]/T)
        if np.isnan(T12):
            T12 = 0
            H12 = 0
    except:
        T12 = 0
        H12 = 0
    try:
        T13 = ascEdges["1_to_3"]/T * np.log2( (ascEdges["1_to_3"]*T)/((ascEdges["1_to_2"]+ascEdges["1_to_3"]+ascEdges["1_to_4"])*(ascEdges["1_to_3"]+ascEdges["2_to_3"])))
        H13 = ascEdges["1_to_3"]/T * np.log2(ascEdges["1_to_3"]/T)
        if np.isnan(T13):
            T13 = 0
            H13 = 0
    except:
        T13 = 0
        H13 = 0
    try:
        T14 = ascEdges["1_to_4"]/T * np.log2( (ascEdges["1_to_4"]*T)/((ascEdges["1_to_2"]+ascEdges["1_to_3"]+ascEdges["1_to_4"])*(ascEdges["1_to_4"] + ascEdges["2_to_4"] + ascEdges["3_to_4"])))
        H14 = ascEdges["1_to_4"]/T * np.log2(ascEdges["1_to_4"]/T)
        if np.isnan(T14):
            T14 = 0
            H14 = 0
    except:
        T14 = 0
        H14 = 0
    try:
        T23 = ascEdges["2_to_3"]/T * np.log2( (ascEdges["2_to_3"]*T)/((ascEdges["2_to_3"]+ascEdges["2_to_4"])*(ascEdges["1_to_3"] + ascEdges["2_to_3"])))
        H23 = ascEdges["2_to_3"]/T * np.log2(ascEdges["2_to_3"]/T)
        if np.isnan(T23):
            T23 = 0
            H23 = 0
    except:
        T23 = 0
        Hc23 = 0
    try:
        T24 = ascEdges["2_to_4"]/T * np.log2( (ascEdges["2_to_4"]*T)/((ascEdges["2_to_4"]+ascEdges["2_to_3"])*(ascEdges["1_to_4"] + ascEdges["2_to_4"] + ascEdges["3_to_4"])))
        H24 = ascEdges["2_to_4"]/T * np.log2(ascEdges["2_to_4"]/T)
        if np.isnan(T24):
            T24 = 0
            H24 = 0
    except:
        T24 = 0
        H24 = 0
    try:
        T34 = ascEdges["3_to_4"]/T * np.log2( (ascEdges["3_to_4"]*T)/(ascEdges["3_to_4"]*(ascEdges["1_to_4"] + ascEdges["2_to_4"] + ascEdges["3_to_4"])))
        H34 = ascEdges["3_to_4"]/T * np.log2(ascEdges["3_to_4"]/T)
        if np.isnan(T34):
            T34 = 0
            H34 = 0
    except:
        T34 = 0
        Hc34 = 0

    A = T12 + T13 + T14 + T23 + T24 + T34
    H = -1*(H12 + H13 + H14 + H23 + H24 + H34)
    #print("T12: ", T12)
    #print("T13: ", T13)
    #print("T14: ", T14)
    #print("T23: ", T23)
    #print("T24: ", T24)
    #print("T34: ", T34)
    #print ("AMI: ", A)
    #print ("Ascendancy: ", A*T)
    #print ("H: ", H)
    
    return A, A*T, T, ascEdges, H, {"1_to_2": T12, "1_to_3": T13, "1_to_4": T14, "2_to_3": T23, "2_to_4": T24, "3_to_4": T34}

def getModuleAssociationCoverage(module):
    thisMod = list(module)
    #methods_df
    #baseNet
    #multNet
    
    thisSubNet = multNet[thisMod].loc[thisMod]
    moduleEdges = msum(thisSubNet)
    ret = {}
    for method in methods_df.keys():
        net = methods_df[method][thisMod].loc[thisMod]
        ret[method] = (int((msum(AND(thisSubNet, net)) / moduleEdges) * 100))/100
    
    directAssoc = (ret["spieceasi"] + ret["ecocopula"]) / 2
    indirectAssoc = (ret["propr"] + ret["ccrepe"] + ret["sparcc"]+ ret["spearman"]) / 4
    return ret, {"direct": directAssoc, "indirect": indirectAssoc}

def getTrLeaveNum(G, j):
    tr = trDict2[j]
    higherNum = 0 # edges to higher trophic level
    lowerNum = 0 # edges from lower trophic level
    for conn in list(G.edges(j)):
        if trDict2[conn[0]] > tr or trDict2[conn[1]] > tr:
            higherNum += 1
        if trDict2[conn[0]] < tr or trDict2[conn[1]] < tr:
            lowerNum += 1
    return lowerNum, higherNum

def getEdgeAMI(i,j, G, T, vabunds):
    #if trDict2[i] > trDict2[j]:
    #    h = i
    #    i = j
    #    j = h
    
    Tij = G[edge[0]][edge[1]]["weight"]
    Ti_out = np.sum([G.get_edge_data(*i_edge)["weight"] for i_edge in G.edges(edge[0])])
    Tj_in = np.sum([G.get_edge_data(*j_edge)["weight"] for j_edge in G.edges(edge[1])])
    if Ti_out == 0 or Tj_in == 0:
        import pdb; pdb.set_trace()
    val = Tij/T * np.log2((Tij*T)/(Ti_out*Tj_in))
    return val

def AMI_taxa(G, vabunds):
    degs = G.degree()
    vabunds = vabunds * 1
    connectedNodes = []
    for edge in list(G.edges):
        connectedNodes.append(edge[0])
        connectedNodes.append(edge[1])
        e0 = vabunds[[edge[0]]].values[0][0] / degs[edge[0]]
        e1 = vabunds[[edge[1]]].values[0][0] / degs[edge[1]]
        G[edge[0]][edge[1]]["weight"] = (e0+e1)/2
    connectedNodes = list(set(connectedNodes))
    T = np.sum(np.sum(vabunds[connectedNodes]))
    
    # import pdb; pdb.set_trace() 
    # T = len(list(G.edges))
    
    AMI = 0
    for edge in list(G.edges):
        #if trDict2[edge[0]] == trDict2[edge[1]]:
        val =  (getEdgeAMI(edge[0],edge[1], G, T, vabunds) + getEdgeAMI(edge[1],edge[0], G, T, vabunds))/2
        #else:
        #    val =  getEdgeAMI(edge[0],edge[1], degs,G,T)
        AMI += val
    print ("AMI: ", AMI)
    print ("Ascendancy: ", AMI*T)
    return AMI, AMI*T

def capacity(G, ascDicts):
    C = 0
    T = len(list(G.edges))
    ascEdges = ascDicts[0]
    for key in ascEdges.keys():
        newC = ascEdges[key] * np.log(ascEdges[key]/T)
        if not np.isnan(newC): 
            C+= ascEdges[key] * np.log(ascEdges[key]/T)
    C = -1 * C
    print("Capacity: ",C)
    return C

def getAscEdges(G):
    ascEdges = {"1_to_2": 0, "1_to_3": 0, "1_to_4": 0, "2_to_3": 0, "2_to_4": 0, "3_to_4": 0}
    nAscEdges = {"1": 0, "2": 0, "3": 0, "4": 0}
    
    ascEdges_list =  {"1_to_2": [], "1_to_3": [], "1_to_4": [], "2_to_3": [], "2_to_4": [], "3_to_4": []}
    for edge in G.edges: 
        le = trDict2[edge[0]]
        re = trDict2[edge[1]]
        if le != re:
            ascEdges[str(min(le,re))+"_to_"+str(max(le,re))] += 1
            ascEdges_list[str(min(le,re))+"_to_"+str(max(le,re))].append(edge)
        if le == re:
            nAscEdges[str(le)] += 1
    return [ascEdges, nAscEdges, ascEdges_list]

def entmischung(net, tops, tneg):
    badPos = []
    badNeg = []
    for k in range(0, len(tpos)):
        if sum(net.loc[tpos[k]])==0:
          badPos.append(tpos[k])
    for k in range(0, len(tneg)):
        if sum(net.loc[tneg[k]])==0:
          badNeg.append(tneg[k])
    
    nTpos = len(tops) - len(badPos)
    nTneg = len(tneg) - len(badNeg)
    
    try:
        tposLinks = msum(net.loc[tpos][tpos])
        nPosLinks = (tposLinks/2)/(((nTpos*nTpos)- nTpos)/2)
    except:
        nPosLinks = 0
    try:
        tnegLinks = (msum(net.loc[tneg][tneg]))
        nNegLinks = (tnegLinks/2)/(((nTneg*nTneg)- nTneg)/2)
    except:
        nNegLinks = 0  
        
    try:
        tmixLinks = msum(net.loc[tpos][tneg])
        nMixedLinks = tmixLinks/(nTneg*nTpos)/2
    except:
        nMixedLinks = 0
    
    print("- pos link rate: ", nPosLinks)
    print("- neg link rate: ", nNegLinks)
    print("- mix link rate: ", nMixedLinks)
    
    return nPosLinks, nNegLinks, nMixedLinks

workdir = 'C:/Users/vdinkel/Desktop/Manuscript/'
metadir = workdir + "output/"
methods = ['propr', 'ccrepe', 'spieceasi', 'esabo', 'ecocopula', 'sparcc', "spearman"]
families_path = workdir+'input/families.csv'

# Read family filter list
df_families = pd.read_table(families_path, sep=";", header=None)
trophicLevels = pd.read_csv(workdir+"input/trophic_levels.csv", delimiter=",", header=None, index_col=None)

# relative abundance data
s_M = pd.read_table(workdir+"output/F_KL-77_rel.csv", sep=";", header=0)

# Temperature correlation files
s_ngrip = pd.read_table(metadir+"KL-77_ngrip.csv", sep=";", header=0)
s_scorrs = pd.read_table(metadir+"KL-77_temp_spearman_corrs.csv", sep=";", header=1)
s_scorrs_p = pd.read_table(metadir+"KL-77_temp_spearman_corrs_p.csv", sep=";", header=1)

# heikes_path = basically workdir
path = workdir+"output/"
s_propr = pd.read_csv(path+"KL-77_propr.csv", delimiter=";", header=0, index_col=0)
s_ccrepe = pd.read_csv(path+"KL-77_ccrepe.csv", delimiter=";", header=0, index_col=0)
s_spieceasi = pd.read_csv(path+"KL-77_spieceasi.csv", delimiter=";", header=0, index_col=0)
s_sparcc = pd.read_csv(path+"KL-77_sparcc.csv", delimiter=";", header=0, index_col=0)
s_esabo = pd.read_csv(path+"KL-77_esabo.csv", delimiter=";", header=0, index_col=0)
#s_ecocopula = pd.read_csv(path+"KL-77_ecocopula_ngrip.csv", delimiter=",", header=0, index_col=0)
s_ecocopula = pd.read_csv(path+"KL-77_ecocopula.csv", delimiter=",", header=0, index_col=0)
s_spearman = pd.read_csv(path+"KL-77_spearman.csv", delimiter=";", header=0, index_col=0)

s_methods_df = {"propr": s_propr, "ccrepe": s_ccrepe, "spieceasi": s_spieceasi,"esabo": s_esabo, "ecocopula": s_ecocopula, "sparcc": s_sparcc, "spearman": s_spearman}

s_sc_df = s_scorrs.T.iloc[1:]
s_sc_pos = list(s_sc_df[s_sc_df > 0.3].dropna().index)
s_sc_neg = list(s_sc_df[s_sc_df < -0.3].dropna().index)

# SUPPLEMENTS: Families with NGRIP correlation
ngrip_corrs = s_sc_df.copy(deep=True)
ngrip_corrs.insert(1, "ngrip_p", list(s_scorrs_p.T.iloc[1:][0]), True)
ngrip_corrs.columns = ["ngrip_r", "ngrip_p"]
ngrip_corrs.to_csv(workdir+"supplementary_information/families_ngrip_corrs.csv", sep=',', index=True, encoding='utf-8')
# SUPPLEMENTS: Families with functions and trophic associations 
fam_func_troph = df_families.copy(deep=True)
fam_func_troph.insert(2, "trophic_level", [int(trophicLevels[trophicLevels[0] == k][1]) for k in list(df_families[1])], True)
fam_func_troph.columns = ["family", "taxonomic_group", "trophic_level"]
fam_func_troph.to_csv(workdir+"supplementary_information/families_trophic_functions.csv", sep=',', index=True, encoding='utf-8')

# Experiment setup
core = "stella"
temptype = "SC"
filename = core+"_"+temptype+"_logic_benchmarks_new.csv"

methods_df = s_methods_df
tpos = s_sc_pos
tneg = s_sc_neg
scorrs = s_scorrs
abunds = s_M
ngrip = s_ngrip

multiplex_A =  methods_df['propr'] + methods_df['spieceasi'] + methods_df['esabo'] + methods_df['sparcc'] + methods_df['spearman'] + methods_df['ecocopula'] + methods_df['ccrepe']
nmult = 7
multiplex_A = multiplex_A / nmult #max(multiplex_A.values.flatten())
combs = {}
weighted_combs = {}
    
method = "spieceasi"
baseNet = methods_df[method]
multNet = AND(methods_df[method], getMultiplexThresh(multiplex_A, 4/7))

mods1 = modcorrs(baseNet, tpos, tneg, 0.0)
mods2 = modcorrs(multNet, tpos, tneg, 0.0)

pre_nPosLinks, pre_nNegLinks, pre_nMixedLinks = entmischung(baseNet, tpos, tneg)
post_nPosLinks, post_nNegLinks, post_nMixedLinks = entmischung(multNet, tpos, tneg)

G_base = nx.from_pandas_adjacency(baseNet)
G_masked = nx.from_pandas_adjacency(multNet)

print("-------------")
print("Before | After")
print("nodes: ", mods1["nodes"], " | ", mods2["nodes"])
print("edges: ", mods1["edges"], " | ", mods2["edges"])
print("modularity: ", mods1["modularity"], " | ", mods2["modularity"])
print("distcomms: ", mods1["distcomms"], " | ", mods2["distcomms"])
print("mixedcomms: ", mods1["mixedcomms"], " | ", mods2["mixedcomms"])
print("--acc/bad--> ",mods1["accoccs"],"/",mods1["badoccs"]," | ",mods2["accoccs"],"/",mods2["badoccs"] )
print("~~~~")
print("posLinkage: ", pre_nPosLinks, " | ", post_nPosLinks, " ~~> ", post_nPosLinks/pre_nPosLinks)
print("negLinkage: ", pre_nNegLinks, " | ", post_nNegLinks, " ~~> ", post_nNegLinks/pre_nNegLinks)
print("mixedLinkage: ", pre_nMixedLinks, " | ", post_nMixedLinks, " ~~> ", post_nMixedLinks/pre_nMixedLinks)

# make cold / warm networks
print ("")
print ("----- Cold | Warm statistics")

thisNet = multNet#baseNet#multNet #baseNet # multNet
thisG = G_masked#G_base#G_masked #G_base # G_masked
thisMods = modcorrs(thisNet, tpos, tneg, 0.0)

trDict = {"1": [], "2": [], "3": [], "4": []} # dictionary with trophic level as key
trDict2 = {} # dictionary with taxon as key and value = trophic level
nodesList = list(thisG.nodes)
for k in nodesList:
    try:
        trGrp = df_families[df_families[0] == k][1].values[0]
        trLevel = trophicLevels[trophicLevels[0] == trGrp][1].values[0]
        trDict[str(trLevel)].append(k)
        trDict2[k] = trLevel
    except:
        import pdb; pdb.set_trace()
    


np.warnings.filterwarnings('ignore')

print ("----")
for method in methods_df.keys():
    print (method, ": ",msum(AND(multNet, methods_df[method])) / msum(multNet))

# HERE OLD PLOT OF NETWORK
import random
def plotNet(ps, ns):
    
    x_offset = 300
    y_offset = 1000
    comms = [list(k) if len(list(k))>1 else [] for k in thisMods['comms']]
    #len(list(set([k[0] for k in thisG.edges])))
    #l0 = (list(set([k[0] for k in thisG.edges])))
    #l1 = (list(set([k[0] for k in thisG.edges])))
    #import pdb; pdb.set_trace()
    plt.figure(figsize=[18,10], dpi = 300)
    pos = nx.spring_layout(thisG, center = [0,0], seed = 6, scale=1500) #4
    #nx.draw_networkx_nodes(G_masked, pos)
    #plt.show()

    trColors = ["#1D8348", "#D4AC0D", "#A04000", "#884EA0"]
    trHeights = 500
    
    i = 0
    flipFlop = 1
    colorBy = "ngrip"#"bio"
    
    '''
    for k in range(0,3):
        if k == 0:
            plt.axhline(y = y_offset*k+trHeights, color = 'r', linestyle = '-')
        else:
            plt.axhline(y = y_offset*k+trHeights, color = 'r', linestyle = '--')
    '''
    plt.axhspan(-500, 500, facecolor='0.2', alpha=0.1)
    #plt.axhspan(500, 1500, facecolor=np.array([176/255,185/255, 58/255]), alpha=0.1)
    plt.axhspan(1500, 2500, facecolor='0.2', alpha=0.1)
    #plt.axhspan(2500, 3500, facecolor=np.array([217/255,33/255, 223/255]), alpha=0.1)
    #plt.axvline(x = 0, color = 'r', linestyle = '-')
    
    nonemptycomms = [com for com in comms if len(com) > 2] #2
    forwardI = [k for k in range(0, len(nonemptycomms))]
    reverseI = [k for k in range(0, len(nonemptycomms))]
    reverseI.reverse()
    if (len(nonemptycomms) % 2) == 0:
        sortedComms = reverseI[::2] + forwardI[::2] #reverseI[::2] + [k+1 for k in forwardI[::2] if k+1 <len(nonemptycomms)]
    else:
        sortedComms = reverseI[::2] + [k+1 for k in forwardI[::2] if k+1 <len(nonemptycomms)]
    
    maxNodeSize = 300
    commPadding = 100
    prevStartX = 0
    prevEndX = 0
    allTRs = {"1": [], "2": [], "3": [], "4": []}
    
    for comm in [nonemptycomms[comI] for comI in sortedComms]: #nonemptycomms
    
        colors = []
        randCol = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])][0]
        trTaxa = {}
        trLevels = {"1": [], "2": [], "3": [], "4": []}
        trGrps = {"1": [], "2": [], "3": [], "4": []}
        
        # get trophic levels of community
        for k in comm:
            trGrp = df_families[df_families[0] == k][1].values[0]
            trLevel = trophicLevels[trophicLevels[0] == trGrp][1].values[0]
            trTaxa[k] = trLevel
            trLevels[str(trLevel)].append(k)
            trGrps[str(trLevel)].append(trGrp)
            allTRs[str(trLevel)].append(k)
        trCounts = Counter([trTaxa[k] for k in trTaxa.keys()])
        maxCounts = trCounts.most_common()[0][1]
        
        node_sizes = dict(np.sum(abunds[comm])*maxNodeSize)
        for key in node_sizes.keys():
            node_sizes[key] = max(min(node_sizes[key], maxNodeSize), maxNodeSize/6)
        moduleSizeSum = sum([node_sizes[k] for k in trLevels[str(trCounts.most_common()[0][0])]])
        
        commStartX = prevEndX
        commEndX = commStartX + sum([node_sizes[k] for k in trLevels[str(trCounts.most_common()[0][0])]])
        commWidth = abs(commEndX - commStartX)
        prevStartX = commStartX
        prevEndX = commEndX + commPadding
        
        ''' X AXIS '''
        # spread across x width of the module
        for tr in trLevels.keys():
            taxaDict = {}
            for tax in trLevels[tr]:
                taxaDict[tax] = pos[tax][0]
            taxaDict = {k: v for k, v in sorted(taxaDict.items(), key=lambda item: item[1])} # sorted
            
            if (len(trLevels[tr])-1) > 0:
                spaceBetween = commWidth / (len(trLevels[tr])-1)
            else:
                spaceBetween = commWidth
            
            for tax in taxaDict.keys():
                try:
                    pos[tax][0] = int(commStartX+(spaceBetween * list(taxaDict.keys()).index(tax)))
                except:
                    import pdb; pdb.set_trace()
        
        ''' Y AXIS '''
        for tr in trLevels.keys():
            y_base =  y_offset*(int(tr)-1)
            y_upper = y_base + trHeights - 100
            y_lower = y_base - trHeights + 100
            
            taxaDict = {}
            for tax in trLevels[tr]:
                taxaDict[tax] = pos[tax][1]
            taxaDict = {k: v for k, v in sorted(taxaDict.items(), key=lambda item: item[1])} # sorted
            
            try:
                spaceBetween = abs(y_upper - y_lower) / (len(trLevels[tr])-1)
            except:
                spaceBetween = trHeights
            for tax in taxaDict.keys():
                pos[tax][1] = int(y_lower+(spaceBetween * list(taxaDict.keys()).index(tax)))
            
            #import pdb; pdb.set_trace()
            
        for k in comm:
            
            trLevel = trTaxa[k]
            '''
            y_base =  y_offset*(trLevel-1)
            y_upper = y_base + trHeights
            y_lower = y_base - trHeights
            
            y_uAttractor = y_upper - (trHeights/1.5)
            y_lAttractor = y_lower + (trHeights/1.5)
            
            #posX = pos[k][0]+(x_offset*i*flipFlop)
            posY = pos[k][1]+y_offset*(trLevel-1)
            
            posY = min(posY, y_upper-40) 
            posY = max(posY, y_lower+40) 
            
            if abs(y_lAttractor - posY) > abs(y_uAttractor - posY): # upper attractor is closer
                y_cAttractor = y_uAttractor
            else:
                y_cAttractor = y_lAttractor
            posY = (posY - y_cAttractor)/1.3 + y_cAttractor
            
            pos[k][1] = posY#np.array([posX, posY])
            '''
            if colorBy == "module":
                colors.append(randCol)  
            if colorBy == "trophic":
                colors.append(trColors[trLevel-1])  
            if colorBy == "ngrip":
                if k in ps:
                    colors.append(np.array([0.9, 0.0, 0.0, 0.85]))  
                elif k in ns:
                    colors.append(np.array([0.0, 0.0, 0.7, 0.85]))
                else:
                    colors.append(np.array([0.7, 0.7, 0.7, 0.7]))
            if colorBy == "bio":
                if k in ps:
                    colors.append(np.array([0.0, 0.39, 0.0, 0.7]))  
                elif k in ns:
                    colors.append(np.array([0.58, 0.0, 0.83, 0.7]))
                else:
                    colors.append(np.array([0.7, 0.7, 0.7, 0.7]))
            
        #plt.axvline(x = x_offset*i*flipFlop, color = 'gray', linestyle = 'dotted')
        
        
        #
        nx.draw_networkx_nodes(thisG, pos, nodelist=comm, node_color=colors, node_size=list(node_sizes.values()))
        
        #nx.draw_networkx_labels(thisG, labels={n: n for n in comm}, pos=pos, font_size=15)
        for node in comm:
            labelSize = max((node_sizes[node] / maxNodeSize)*13, 7)
            plt.text(pos[node][0], pos[node][1], node, fontsize=labelSize, ha='center', va='center')
        i+=1
        flipFlop *= -1
        
        #import pdb; pdb.set_trace()
    
    edgecols = []
    alphas = []
    for edge in thisG.edges:
        baseCol = np.array([0.7, 0.7, 0.7, 0.7])
        if colorBy == "ngrip":
            if edge[0] in ps or edge[1] in ps:
                baseCol = np.array([0.9, 0.0, 0.0, 0.7])
            if edge[0] in ns or edge[1] in ns:
                baseCol = np.array([0.0, 0.0, 0.9, 0.7])
            if (edge[0] in ps and edge[1] in ns) or (edge[0] in ns and edge[1] in ps):
                baseCol = np.array([60/255, 0.0, 100/255, 0.7])
        if colorBy == "bio":
            if edge[0] in ps or edge[1] in ps:
                baseCol = np.array([0.0, 0.9, 0.0, 0.85])
            if edge[0] in ns or edge[1] in ns:
                baseCol = np.array(np.array([0.58, 0.0, 0.83, 0.85]))
            if (edge[0] in ps and edge[1] in ns) or (edge[0] in ns and edge[1] in ps):
                baseCol = np.array([1.0, 0.0, 0.65, 1.0])
        edgecols.append(baseCol)
        alphas.append(multiplex_A[edge[0]][edge[1]])
        #import pdb; pdb.set_trace()
    
    nx.draw_networkx_edges(thisG, pos, edge_color = edgecols, width=np.array(alphas)*2, alpha=np.array(alphas))
    plt.xlim(0-150, commEndX+150)
    plt.ylim(-500, 3500)
    import time
    obj = time.gmtime(0)
    curr_time = round(time.time()*1000)
    plt.savefig("C:/Users/vdinkel/Desktop/PhD/Presentations/networkVisialization/"+str(curr_time)+".png")
    plt.show()
    plt.close()
#plotNet()

#import pdb; pdb.set_trace()
#for i in range(0,len(abunds)):
#plotNetwork(False)
G_new = nx.from_pandas_adjacency(multNet)
degrees = list(G_new.degree)
for deg in degrees:
    if deg[1] == 0:
        G_new.remove_node(deg[0])
        #print("removed: ", deg[0])
exportToGephi(G_new, df_families, tpos, tneg, thisMods, path+"cn_spieceasi_05.gml", trophicLevels)

# TROPHIC STATISTICS
# 1 = 17.1; 2 = 3.5; 3 = 20.5; 4 = 0.9
abunds[[k for k in trDict2.keys() if trDict2[k] == 1]].sum().sum()


# COMMUNITY STATISTICS
lcc = list(max(nx.connected_components(G_new), key=len))
between_keystones = nx.betweenness_centrality(G_new, k=None, normalized=True, weight=None)
degree_keystones = nx.eigenvector_centrality(G_new)

df_modules = {}
df_keystones = {}

i = 1
for thiscom in thisMods["comms"]:
    thiscom = list(thiscom)
    if thiscom[0] in lcc:
        mtl, rmtl = getModuleMTL(thiscom)
        tempcorr = getModuleTemp(thiscom)
        df_modules[i] = {"size": len(thiscom), "functional_composition": df_families[df_families[0].isin(thiscom)][1].value_counts(), "mean_trophic_level": mtl, "interglacial": tempcorr['pos'], "glacial": tempcorr['neg']}
        for fam in thiscom:
            df_keystones[fam] = {"module": i, "degree_centrality": degree_keystones[fam], "betweenness_centrality": between_keystones[fam], "keystonness": (degree_keystones[fam] +  between_keystones[fam])}
        i += 1
df_modules = pd.DataFrame(df_modules)
df_keystones = pd.DataFrame(df_keystones).T

df_modules.to_csv(workdir+"supplementary_information/cn_modules.csv", sep=';', index=True, encoding='utf-8')
df_keystones.to_csv(workdir+"supplementary_information/cn_centralities.csv", sep=';', index=True, encoding='utf-8')

'''
import plotly.graph_objects as go
gephi_export_file = "C:/Users/vdinkel/Desktop/Manuscript/output/gephi_data_export.csv"
gephi_data = pd.read_table(gephi_export_file, sep=";")
community_file = "C:/Users/vdinkel/Desktop/Data/Manuscript/all_communities.csv"
community_data = pd.read_table(community_file, sep=";")
rank = "order"
for comm_i in range(0,4):
    print("Comm: ", comm_i)
    print(community_data[community_data["community"] == comm_i][rank].value_counts())
import copy

import plotly.graph_objects as go
import matplotlib.colors as pltclr
# = ["0", "0", "1", "1", "0"]
#target = [2, 3, 4, 5, 4]
#value = [8, 2, 2, 8, 4]

trColors = ["#009988", "#8ad963", "#0077BB", "#EE3377"] #33BBEE
colors = {
    "primary_producer": pltclr.to_rgba(trColors[0]),
    "primary_consumer": pltclr.to_rgba(trColors[1]),
    "secondary_consumer": pltclr.to_rgba(trColors[2]),
    "tertiary_consumer": pltclr.to_rgba(trColors[3]),
    "module1": pltclr.to_rgba("#5CD3FF"), #0077BB
    "module234": pltclr.to_rgba("#FF6644"), #FF3F14 ##CC3311
    }

for key in colors.keys():
    colors[key] = tuple([int(colors[key][0] * 255), int(colors[key][1] * 255), int(colors[key][2] * 255), 1.0])
test = list(colors.keys())
for key in test:
    colors[key+"_l"] =tuple([colors[key][0], colors[key][1], colors[key][2], 0.6])

for key in colors.keys():
    colors[key] = 'rgba'+str(colors[key])

label_id = {}
id_label = {}

id_count = 8

id_label[4] = "1: Primary Producer"
label_id["1: Primary Producer "] = 4
id_label[5] = "2: Primary Consumer"
label_id["2: Primary Consumer"] = 5
id_label[6] = "3: Secondary Consumer"
label_id["3: Secondary Consumer"] = 6
id_label[7] = "4: Tertiary Consumer"
label_id["4: Tertiary Consumer"] = 7

for module_i in range(0,4):
    id_label[module_i] = "Module "+str(module_i)
    label_id["Module "+str(module_i)] = module_i
    
    value_counts = gephi_data[gephi_data["module"]==module_i]["taxonomic_group"].value_counts()
    for k in range(0, len(value_counts.index)):
        if value_counts.index[k] not in label_id.keys():
            label_id[value_counts.index[k]] = id_count
            id_label[id_count] = value_counts.index[k]
            id_count += 1

source = []
target = []
value = []
link_cols = []
trophic_levels_f = "C:/Users/vdinkel/Desktop/Data/Manuscript/input/trophic_levels.csv"
trophic_levels = pd.read_table(trophic_levels_f, sep=",", header = None)
# MODULE <-> TAXONOMIC GROUP
for module_i in range(0,4):
    value_counts = gephi_data[gephi_data["module"]==module_i]["taxonomic_group"].value_counts()
    for i in range(0, len(value_counts)):
        source.append(module_i)
        target.append(label_id[value_counts.index[i]])
        value.append(value_counts[i])
        
        troph_target = trophic_levels[trophic_levels[0] == value_counts.index[i]][1].values[0]
        if troph_target == 1:
            link_cols.append(colors["primary_producer_l"])
        if troph_target == 2:
            link_cols.append(colors["primary_consumer_l"])
        if troph_target == 3:
            link_cols.append(colors["secondary_consumer_l"])
        if troph_target == 4:
            link_cols.append(colors["tertiary_consumer_l"])

# TAXONOMIC GROUP <-> TROPHIC LEVEL
id_taxstart = id_count+1
taxgrp_counts = dict(gephi_data["taxonomic_group"].value_counts())

for i in range(0,len(trophic_levels)):
    try:
        source.append(label_id[trophic_levels[0][i]])
        target.append(trophic_levels[1][i]+3)
        value.append(taxgrp_counts[trophic_levels[0][i]])
        
        if trophic_levels[1][i]+3 == 4:
            link_cols.append(colors["primary_producer_l"])
        if trophic_levels[1][i]+3 == 5:
            link_cols.append(colors["primary_consumer_l"])
        if trophic_levels[1][i]+3 == 6:
            link_cols.append(colors["secondary_consumer_l"])
        if trophic_levels[1][i]+3 == 7:
            link_cols.append(colors["tertiary_consumer_l"])
            
        #import pdb; pdb.set_trace()
    except:
        print ("catched: ", trophic_levels[0][i])

node_colors = [colors["module1"], colors["module234"], colors["module234"], colors["module234"], colors["primary_producer"], colors["primary_consumer"], colors["secondary_consumer"], colors["tertiary_consumer"]]

for lab in [id_label[k] for k in sorted(id_label)][8:]:
    trlvl = trophic_levels[trophic_levels[0] == lab][1].values[0]
    if trlvl == 1:
        node_colors.append(colors["primary_producer_l"])
    if trlvl == 2:
        node_colors.append(colors["primary_consumer_l"])
    if trlvl == 3:
        node_colors.append(colors["secondary_consumer_l"])
    if trlvl == 4:
        node_colors.append(colors["tertiary_consumer_l"])


# SANKEY DIAGRAM EXPORT
fig = go.Figure(data=[go.Sankey(
    orientation = "v",
    node = dict(
      pad = 20,
      thickness = 80,
      #y= [0.5 for k in range(0,len(source))],      #x= [0.2 for i in range(0,4)]+[0.2 for i in range(0,4)]+[0.5 for k in range(8,len(source))],
      line = dict(color = "black", width = 1.1),
      #label = [id_label[k] for k in sorted(id_label)],
      color = node_colors
    ),
    link = dict(
      source = source, # indices correspond to labels, eg A1, A2, A1, B1, ...
      target = target,
      value = value,
      color = link_cols
  ))])

fig.update_layout(
    hovermode = 'x',
    title="<b>Module Composition</b><br />Distribution of Functional Groups and Trophic Levels in Major LCC Modules",
    font=dict(size = 12, color = 'black'),
    plot_bgcolor='white',
    paper_bgcolor='white',
    )

fig.write_html("C:/Users/vdinkel/Desktop/Data/Manuscript/output/sankey-diagram-plotly_all.html")
'''

# ROBUSTNESS
from random import sample
import copy
import operator

def getTrophicLevelDiversity(G):
    retlen = 0
    try:
        retlen = len(list(set(fam_func_troph[fam_func_troph["family"].isin(list(max(nx.connected_components(G), key=len)))]["trophic_level"])))
    except:
        retlen = 0
    return retlen

from collections import Counter
import math

def getRobustnessMetric(G):
    try:
        metric = len(G.nodes)
        #metric = nx.node_connectivity(G)
        #metric = getTrophicLevelDiversity(G)
        #metric = calculate_pielou_evenness(G)
    except:
        import pdb; pdb.set_trace()
        metric = 0
    return metric

#from networkx.algorithms.flow import shortest_augmenting_path
def getMinimumCut(G):
    delVals = {}
    init_lcc = len(list(max(nx.connected_components(G), key=len)))
    for node in G.nodes:
        tempG =  G.copy()
        tempG.remove_node(node)
        delVals[node] = len(list(max(nx.connected_components(tempG), key=len)))
    minVal = min(delVals.values())
    keys = [k for k, v in delVals.items() if v == minVal]
    return keys
        

def minCut(G, normalized, critpath = []):
    minCuts = []
    init_lcc = list(max(nx.connected_components(G), key=len))
    tempG =  G.subgraph(init_lcc).copy()
    init_metric = getRobustnessMetric(tempG)
    if normalized:
        minCuts.append(('', init_metric/init_metric, 'black'))
    else:
        minCuts.append(('', init_metric, 'black'))
    
    for i in range(0, len(G.nodes)):
        try:
        #min_node_cut = list(nx.minimum_node_cut(tempG, flow_func=None))
            if len(critpath) == 0:
                min_node_cut = getMinimumCut(tempG)
            else:
                min_node_cut = [critpath[i]]
            if len(min_node_cut)>1:
                # 1. prefer the removal of highest degree node
                degs = {}
                for cut in min_node_cut:
                    degs[cut] = tempG.degree[cut]
                min_node_cut = [k for k, v in degs.items() if v == max(degs.values())]
                
                # 2. prefer the removal of highest trophic node
                if len(min_node_cut)>1:
                    trophs = {}
                    for cut in min_node_cut:
                        trophs[cut] = trDict2[cut]
                    min_node_cut = [k for k, v in trophs.items() if v == max(trophs.values())]
                
                if len(min_node_cut)>1:
                    min_node_cut = min_node_cut[0]
                if len(min_node_cut)==1:
                    min_node_cut = min_node_cut[0]
                
            else:
                min_node_cut = min_node_cut[0]
            if min_node_cut in tpos:
                col = "#CC3311"
            elif min_node_cut in tneg:
                col = "#0077BB"
            else:
                col = "black"
            tempG.remove_node(min_node_cut)
            lcc = list(max(nx.connected_components(tempG), key=len))
            tempG =  tempG.subgraph(lcc).copy()
            metric = getRobustnessMetric(tempG)
            if normalized:
                metric = metric / init_metric
            minCuts.append((min_node_cut, metric, col))
        except:
            return minCuts
    return minCuts

def getRobustness(G, module, ndels, preferred = False, preflist = [], normalized = False):
    all_metrics = []
    for i in range(0, ndels):
        # Initialize Robustness Parameters
        tempG = G.subgraph(module).copy()
        init_metric = getRobustnessMetric(tempG)
        
        if not normalized:
            metrics = [init_metric]
        else:
            metrics = [init_metric/init_metric]
        if preferred:
            tmp_preflist = preflist.copy()
            
        # Iterative Node Extinction
        for k in range(0, len(G.nodes)):
            try:
                # Decide which node to remove & remove node
                if preferred:
                    del_prefs = []
                    for tmp in tmp_preflist:
                        if tmp not in tempG.nodes():
                            del_prefs.append(tmp)
                    for del_pref in del_prefs:
                        tmp_preflist.remove(del_pref)
                    if len(tmp_preflist)>0:
                        node = sample(tmp_preflist, 1)[0]
                        tmp_preflist.remove(node)
                        if node in list(tempG.nodes):
                            tempG.remove_node(node)
                    else:
                        break
                else:
                    node = sample(list(tempG.nodes()), 1)[0]
                    tempG.remove_node(node)
                
                # get new lcc & get lcc metric
                tempG =  tempG.subgraph(list(max(nx.connected_components(tempG), key=len))).copy()
                metric = getRobustnessMetric(tempG)
                if not normalized:
                    metrics.append(metric) 
                else:
                    metrics.append(metric/init_metric) 
            except:
                metrics.append(0)
                break
        all_metrics.append(metrics)
    
    means = []
    stds = []
    pd_all_metrics = pd.DataFrame(all_metrics)#.fillna(method='ffill', inplace=True)#.fillna(0)
    pd_all_metrics = pd_all_metrics.T.fillna(method='ffill').T
    for k in pd_all_metrics.columns:
        k_vals = pd_all_metrics[k].values
        means.append(np.mean(k_vals))
        stds.append(np.std(k_vals))
    #import pdb; pdb.set_trace()
    return np.array(means), np.array(stds)
            
def moduleRobustness(G, module, ndels, normalized = False, critpath = []):
    
    thisPos = list(set(module).intersection(tpos))
    thisNeg = list(set(module).intersection(tneg))
    thisTemps = thisPos + thisNeg
    
    #lcc = list(max(nx.connected_components(G), key=len))
    rand_means, rand_stds = getRobustness(G, module, ndels, normalized = normalized)
    temp_means, temp_stds = getRobustness(G, module, ndels, True, thisTemps, normalized = normalized)
    minCuts = minCut(G.subgraph(module).copy(), normalized = normalized, critpath = critpath)
    
    #import pdb; pdb.set_trace()
    
    if normalized:    
        xrands = [x/len(module) for x in range(0, len(rand_means))]
        xtemps = [x/len(module) for x in range(0, len(temp_means))]
        xcuts = [x/len(module) for x in range(0, len(minCuts))]
    else:
        xrands = [x for x in range(0, len(rand_means))]
        xtemps = [x for x in range(0, len(temp_means))]
        xcuts = [x for x in range(0, len(minCuts))]
    #import pdb; pdb.set_trace()
    plt.figure(figsize = (7, 5))
    
    import scipy.stats as ss
    dd = {'x': xcuts, 'y': [k[1] for k in minCuts]}
    df_dd = pd.DataFrame(data=dd)
    lm = smf.ols(formula='y ~ x', data=df_dd).fit()
    xlm = [k/10 for k in range(0,10)]
    ylm = [k * lm.params["x"] + lm.params["Intercept"] for k in xlm]
    print("LM: ", lm.params)
    #import pdb; pdb.set_trace()
    plt.fill_between(xrands, rand_means + rand_stds, rand_means - rand_stds, alpha=0.1, color = "gray")
    #plt.fill_between(xtemps, temp_means + temp_stds, temp_means - temp_stds, alpha=0.1, color = "red")
    plt.plot(xrands, rand_means, linestyle="-", label="random targeted", marker = 'o', alpha=0.4, color="gray")
    plt.plot(xlm, ylm, label="critical path slope", alpha=0.2, linewidth= 5, color="black")
    plt.plot([0,1],[1,0], linestyle="--", label="maximum robustness", alpha=0.5, color="gray")
    #plt.plot(xtemps, temp_means, linestyle="-", label="NGRIP targeted", color="red")
    plt.plot(xcuts, [k[1] for k in minCuts], label="critical path", marker = "X", color="black")
    
    i = 0
    for label, value, color in minCuts:
        #t = Affine2D().scale(skill).rotate_deg(takeoff)
        #m = MarkerStyle(SUCCESS_SYMBOLS[mood], transform=t)
        if i == 0:
            m = "*"
            s = 10
        else:
            if color == "#CC3311" or color == "#0077BB" :
                m = "X"
                s = 10
            else:
                m = "X"
                s = 10
        plt.text(xcuts[i] + 0.02, value, label, color = color, fontsize=12)
        plt.plot(xcuts[i], value, marker=m, color=color, markersize = s)
        i += 1
    plt.xlabel("family knockouts (%)")
    plt.ylabel("families remaining (%)")
    if normalized:
        plt.xlim(0, 1)
        plt.ylim(0, 1.01)
    plt.legend()
    # supplements: robustness mincuts
    df_mincuts = pd.DataFrame(minCuts)
    df_mincuts.columns = ["extinction", "families_remaining", "color"]
    df_mincuts.index = xcuts
    df_mincuts.index.name = "knockouts"
    df_mincuts.to_csv(workdir+"supplementary_information/mincuts_module_"+str(delmod+1)+".csv", sep=';', index=True, encoding='utf-8')
    plt.savefig(workdir+'plots/robustness_module_'+str(delmod+1)+'.png', format='png', dpi=1200)
    plt.show()
    
    return rand_means, rand_stds

def getLCCRobustness(G, module, ndels, normalized = False):
    rand_means, rand_stds = getRobustness(G, module, ndels, normalized = normalized)
    #temp_means, temp_stds = getRobustness(G, module, ndels, True, thisTemps, normalized = normalized)
    minCuts = minCut(G.subgraph(module).copy(), normalized = normalized)
    if normalized:    
        xrands = [x/len(module) for x in range(0, len(rand_means))]
        #xtemps = [x/len(module) for x in range(0, len(temp_means))]
        xcuts = [x/len(module) for x in range(0, len(minCuts))]
    else:
        xrands = [x for x in range(0, len(rand_means))]
        #xtemps = [x for x in range(0, len(temp_means))]
        xcuts = [x for x in range(0, len(minCuts))]
    
    
    plt.figure(figsize = (7, 5))
    plt.fill_between(xrands, rand_means + rand_stds, rand_means - rand_stds, alpha=0.1, color = "gray")
    #plt.fill_between(xtemps, temp_means + temp_stds, temp_means - temp_stds, alpha=0.1, color = "red")
    plt.plot(xrands, rand_means, linestyle="-", label="random targeted", marker = 'o', alpha=0.4, color="gray")
    plt.plot([0,1],[1,0], linestyle="--", label="maximum robustness", alpha=0.5, color="gray")
    #plt.plot(xtemps, temp_means, linestyle="-", label="NGRIP targeted", color="red")
    plt.plot(xcuts, [k[1] for k in minCuts], label="critical path", marker = "X", color="black")
    
    i = 0
    for label, value, color in minCuts:
        if i == 0:
            m = "*"
            s = 10
        else:
            if color == "#CC3311" or color == "#0077BB" :
                m = "X"
                s = 10
                plt.text(xcuts[i] + 0.02, value, label+" (module "+str(int(df_keystones.loc[label]["module"]))+")", color = color, fontsize=12)
                #import pdb; pdb.set_trace()
            else:
                m = "X"
                s = 10
                #if (normalized and (minCuts[i-1][1] - value) > 0.05) or ((minCuts[i-1][1] - value) > 5):
                #    plt.text(xcuts[i] + 0.02, value, label+" (module "+str(int(df_keystones.loc[label]["module"]))+")", color = color, fontsize=12)
        
        plt.plot(xcuts[i], value, marker=m, color=color, markersize = s)
        i += 1
    plt.xlabel("family knockouts (%)")
    plt.ylabel("families remaining (%)")
    if normalized:
        plt.xlim(0, 1)
        plt.ylim(0, 1.01)
    plt.legend()
    # supplements: robustness mincuts
    df_mincuts = pd.DataFrame(minCuts)
    df_mincuts.columns = ["extinction", "families_remaining", "color"]
    df_mincuts.index = xcuts
    df_mincuts.index.name = "knockouts"
    df_mincuts.to_csv(workdir+"supplementary_information/mincuts_lcc.csv", sep=';', index=True, encoding='utf-8')
    plt.savefig(workdir+'plots/robustness_lcc.png', format='png', dpi=1200)
    plt.show()
    return rand_means, rand_stds

currmean = []
class NonBinTree:
    def __init__(self, label, G, depth, parent):
        self.label = label
        self.G = G.copy()
        self.nodes = []
        self.lcc = []
        self.leaf = False
        self.depth = depth
        self.parent = parent
        
        if self.label in G.nodes():
            self.rem_node(self.label)
        
        self.lcc = self.get_lcc()
        self.G = self.G.subgraph(self.lcc).copy()
        
        if len(self.lcc) > 0:
            if self.depth > 3 and len(self.lcc) > 7:
                return
            for node in self.lcc:
                self.a_node(node, self.G, self.depth+1, self)
        else:
            self.leaf = True
        
    def rem_node(self, node):
        self.G.remove_node(node)

    def get_lcc(self):
        try:
            lcc = list(max(nx.connected_components(self.G), key=len))
        except:
            lcc = []
        return lcc

    def a_node(self, label, G, depth, parent):
        self.nodes.append(NonBinTree(label, G, depth, parent))
        
    def printPath(self):
        print ("[",self.depth,"] ", self.label," ", len(self.lcc))
        self.parent.printPath()
    
    def getLeafAtDepth(self, depth):
        if self.leaf and self.depth==depth:
            print (self.label," ", self.depth)
            these_leafs.append(self)
            return self

        for node in self.nodes:
            node.getLeafAtDepth(depth)
        #return retlist
    
    def getMinAtDepth(self, depth):
        global currmean
        
        if self.depth == depth and len(currmean) == 0:
            currmean.append(self)
        if self.depth == depth and len(currmean) > 0 and len(currmean[-1].lcc) >= len(self.lcc):
            currmean.append(self)
        
        for node in self.nodes:
            node.getMinAtDepth(depth)
        

    #def __repr__(self):
    #    return f"NonBinTree({self.lcc}): {self.nodes}"
import itertools
#testCom = list(thisMods['comms'][5])
#robG = nx.subgraph(thisG, testCom).copy()
delpath_0 = ["Blenniidae", "Monodontidae", "Salmonidae", "Lateolabracidae", "Batrachoididae", "Petromyzontidae", "Sparidae"]
delpath_2 = ["Phaeocystaceae", "Attheyaceae", "Bacillariaceae", "Ulnariaceae", "Eunotiaceae", "Gomphonemataceae", "Oxystominidae"]
delpath_5 = ["Clupeidae", "Balaenopteridae", "Mytilidae", "Delphinidae"]
#for dpath in delpath:
#    robG.remove_node(dpath)
#    this_lcc = list(max(nx.connected_components(robG), key=len))
#    print(dpath, len(this_lcc))
#    robG = nx.subgraph(robG, this_lcc).copy()

#tree = NonBinTree("", robG, 0, None)
#these_leafs = []
#tree.getLeafAtDepth(4)
#tree.getMinAtDepth(3)
#import pdb; pdb.set_trace()

runRobustness = False
i = 0
if runRobustness:
    del_modules = [0, 2, 5]
    del_critpaths = [delpath_0, delpath_2, delpath_5]
    nruns = 1000
    for delmod in del_modules:
        testCom = list(thisMods['comms'][delmod])
        lcc_G = nx.subgraph(thisG, list(max(nx.connected_components(thisG), key=len))).copy()
        rand_means, rand_stds = moduleRobustness(lcc_G, testCom, nruns, normalized = True, critpath = del_critpaths[i])
        # rand_means, rand_stds = moduleRobustness(thisG, list(thisMods['comms'][0]), 100)
        i += 1
    lcc_G = nx.subgraph(thisG, list(max(nx.connected_components(thisG), key=len))).copy()
    getLCCRobustness(lcc_G, lcc_G.nodes, nruns, normalized = True)
#import pdb; pdb.set_trace()
AMIs = []
AMI_comps = []
As= []
A_comps = []
Ts = []
ascEdgesL = []
Hs = []
T_states_L = []

def getFamMod(fam, mods):
    for i in range(0,7):
        if fam in list(mods[i]):
            return i
    return None

def calcAsc(G, trDict2, whichAscs):
    ascs = {}
    ascs_mods = {}
    for i in range(0,len(mods2["comms"])):
        ascs_mods[str(i)] = 0.0
    T =sum([k[2]["weight"] for k in G.edges(data=True)])
    AMI = 0.0
    Hc = 0.0
    H = 0.0
    DC = 0.0

    all_ascs = {}
    all_ins = {}
    all_outs = {}
    for ff in range(1,5):
        for fff in range(1,5):
            ascs[str(ff)+"_to_"+str(fff)] = 0.0
            all_ascs[str(ff)+"_to_"+str(fff)] = []
    
    for edge in G.edges(data=True):
        Tj = 0.0
        Ti = 0.0
        for inedge in G.in_edges(edge[1], data=True):
            Tj += inedge[2]["weight"]
        for outnedge in G.out_edges(edge[0], data=True):
            Ti += outnedge[2]["weight"]   
        Tij = edge[2]["weight"]
        
        thisAMI = Tij/T * np.log2((Tij * T)/(Ti * Tj))
        thisA = thisAMI * T
        thisHc = -1 * (Tij/T * np.log2((Tij*Tij)/T))
        thisH = thisAMI + thisHc
        thisDC = T * thisH
        
        if whichAscs == "A/DC":
            thisADC =  thisA / thisDC 
        if whichAscs == "A":
            thisADC =  thisA
        if whichAscs == "DC":
            thisADC =  thisDC
        if whichAscs == "H":
            thisADC =  thisH
        
        AMI += thisAMI
        Hc += thisHc
        H += thisH
        DC += thisDC
        
        if edge[0] not in all_outs.keys():
            all_outs[edge[0]] = []
        if edge[1] not in all_ins.keys():
            all_ins[edge[1]] = []
        all_outs[edge[0]].append(thisADC)
        all_ins[edge[1]].append(thisADC)

        if str(trDict2[edge[0]])+"_to_"+str(trDict2[edge[1]]) not in ascs.keys():
            ascs[str(trDict2[edge[0]])+"_to_"+str(trDict2[edge[1]])] = thisADC
            all_ascs[str(trDict2[edge[0]])+"_to_"+str(trDict2[edge[1]])] = [thisADC]
        else:
            ascs[str(trDict2[edge[0]])+"_to_"+str(trDict2[edge[1]])] +=  thisADC
            all_ascs[str(trDict2[edge[0]])+"_to_"+str(trDict2[edge[1]])].append(thisADC)
        
        e0_mod = getFamMod(edge[0], mods2["comms"])
        e1_mod = getFamMod(edge[1], mods2["comms"])

        try:
            if (e0_mod == e1_mod) and (e0_mod != None):
                ascs_mods[str(e0_mod)] += thisADC
        except:
            import pdb; pdb.set_trace()
            
        #import pdb; pdb.set_trace()

    return T, AMI, Hc, H, ascs, all_ins, all_outs, ascs_mods

def getSimG(G, vabunds, trDict2):
    diG = nx.DiGraph()
    # NORMALIZE ABUNDANCE VECTOR TO 1 ? (this sets system throughput to constant)
    #if vabunds[list(G.nodes)].sum().sum() != 1.0:
    #    vabunds = vabunds[list(G.nodes)] / vabunds[list(G.nodes)].sum().sum()
    nodes = list(G.nodes())
    for node in nodes:
        if vabunds[node].values[0]>0.0:
            diG.add_node(node, level = trDict2[node])
    diNodes = list(diG.nodes())
    edges = list(G.edges())
    for edge in edges:
        if edge[0] in diNodes and edge[1] in diNodes:
            # transform an edge from G into two directed edges in diG
            diG.add_edge(edge[0], edge[1])
            diG.add_edge(edge[1], edge[0])
    degs = diG.out_degree()
    for edge in list(diG.edges()):
        # loop through all directed edges in diG, divide abundance by out degree & set edge weight
        diG[edge[0]][edge[1]]["weight"] = vabunds[edge[0]].values[0] / degs[edge[0]]
    return diG

'''
# NEW ASCENDENCY WITH DIRECTED NETWORK
'''
def flowIndeces(abunds, proxievals, trDict2, flowG, whichAscs):
    # proxievals = ngrip["temp"].values
    AMIs = []
    Ts = []
    Hs = []
    ascss = {}
    all_flows = []
    all_mods_flows = []
    lccs = {"size": [], "edges": []}
    for k in range(0,len(abunds)):
        vabunds = abunds.iloc[[k]]
        simG = getSimG(flowG, vabunds, trDict2)
        
        T, AMI, Hc, H, ascs, all_ins, all_outs, ascs_mods =  calcAsc(simG, trDict2, whichAscs)
        AMIs.append(AMI)
        Ts.append(T)
        Hs.append(H)
        all_flows.append([all_outs, all_ins])
        
        for key in ascs.keys():
            if key not in ascss:
                ascss[key] = [ascs[key]]
            else:
                ascss[key].append(ascs[key])
        
        all_mods_flows.append(ascs_mods)
        lccs["size"].append(len(max(nx.connected_components(nx.DiGraph.to_undirected(simG)), key=len)))
        lccs["edges"].append(len(nx.DiGraph.to_undirected(simG).edges()))
    
    As = np.array(AMIs) * np.array(Ts)
    DCs = np.array(Ts) * np.array(Hs)
    print("-----------")
    print("AMI:  -> ", scipy.stats.pearsonr(AMIs, proxievals))
    print("T:  -> ", scipy.stats.pearsonr(Ts, proxievals))
    print("A: -> ", scipy.stats.pearsonr(As, proxievals))
    print("H: -> ", scipy.stats.pearsonr(Hs, proxievals))
    print("DC:  -> ", scipy.stats.pearsonr(DCs, proxievals))
    print("A/DC:  -> ", scipy.stats.pearsonr(As/DCs, proxievals))
    
    maturity = {"AMI": AMIs, "T": Ts, "A": As, "H": Hs, "DC" : DCs, "A/DC": As/DCs }
    retAscss = {}
    
    for key in ascss:
        retAscss[key] = {"r": 0.0, "p": 0.0}
        r,p = scipy.stats.pearsonr(ascss[key], proxievals)
        if p <= 0.05:
            print (key, " -> ", r, " p = ",p)
            retAscss[key]["r"] = r
            retAscss[key]["p"] = p
    
    return retAscss, ascss, all_flows, all_mods_flows, lccs, maturity

# SPIECEASI WEIGHTED TEST
performHTComparison = False
if performHTComparison:
    s_spieceasi_weighted = pd.read_csv(path+"KL-77_spieceasi_weighted.csv", delimiter=";", header=0, index_col=0)
    s_spieceasi_edgelist = pd.read_csv(path+"KL-77_spieceasi_weighted_edgelist.csv", delimiter=";", header=0, index_col=0)
    #se_weighted = getMultiplexThresh(s_spieceasi_weighted, 0.3421) #0.327
    #mods3 = modcorrs(se_weighted, tpos, tneg, 0.0, globi_M)
    
    ht_net = s_spieceasi_weighted
    
    mods_base = []
    mods_ht = []
    
    entmisch_base = []
    entmisch_ht = []
    
    alreadyPassed = []
    all_compres = []
    compres = {}
    for i in range(0, 9):
        print (i)
        j = i/10
        #j = 0.57
        cn_j = AND(baseNet, getMultiplexThresh(multiplex_A, j))
    
        print("--mult: ", str(msum(cn_j)))
        if msum(cn_j) not in alreadyPassed:
            modsMult = modcorrs(cn_j, tpos, tneg, 0.0)
            mods_base.append(modsMult)
            
            enmisch_b = entmischung(cn_j, tpos, tneg)
            entmisch_base.append(enmisch_b)
            
            tmp_net = []
            #prevPass = 0
            for k in range(0, 1000):
                if msum(getMultiplexThresh(ht_net, k/1000)) <= msum(cn_j):
                    tmp_net = getMultiplexThresh(ht_net, k/1000)
                    break
                #else:
                #    prevPass = k
        
            # cn_j.sum()[cn_j.sum() > 0.0] # amount of nodes with degree > 0
            print("--ht: ", str(msum(tmp_net)))
            
            modsHT = modcorrs(tmp_net, tpos, tneg, 0.0)
            mods_ht.append(modsHT)
            entmisch_h = entmischung(tmp_net, tpos, tneg)
            entmisch_ht.append(entmisch_h)
            
            G1 = nx.from_pandas_adjacency(cn_j)
            G2 = nx.from_pandas_adjacency(tmp_net)
            
            g1_lcc = len(max(nx.connected_components(G1), key=len))
            g2_lcc = len(max(nx.connected_components(G2), key=len))
            
            print("G1 edges: ", len(G1.edges))
            print("G2 edges: ", len(G2.edges))
        
            #flow_G1, flow_G1_ascss, all_flows_G1 = flowIndeces(abunds, ngrip, trDict2, G1)
            #flow_G2, flow_G2_ascss, all_flows_G2  = flowIndeces(abunds, ngrip, trDict2, G2)
        
            #g1_comps = [len(c) for c in sorted(nx.connected_components(G1), key=len, reverse=True)]
            #g2_comps = [len(c) for c in sorted(nx.connected_components(G2), key=len, reverse=True)]
            
            alreadyPassed.append(msum(cn_j))
    
            compres[str(j)] = {
                       "cn_edges" : len(G1.edges), "cn_modularity" : modsMult["modularity"], "cn_cc" : modsMult["cc"], "cn_lcc": g1_lcc,"cn_distinct_comms" : modsMult["distcomms"], "cn_mixed_comms" : modsMult["mixedcomms"], "cn_pos_linkage" : enmisch_b[0], "cn_neg_linkage" : enmisch_b[1], "cn_mix_linkage" : enmisch_b[2],
                       "ht_edges" : len(G2.edges), "ht_modularity" : modsHT["modularity"], "ht_cc" : modsHT["cc"], "ht_lcc": g2_lcc,"ht_distinct_comms" : modsHT["distcomms"], "ht_mixed_comms" : modsHT["mixedcomms"], "ht_pos_linkage" : entmisch_h[0], "ht_neg_linkage" : entmisch_h[1], "ht_mix_linkage" : entmisch_h[2]}
            #all_compres.append(compres)
        else:
            pass
    # SUPPLEMENTS: comparison of base network (cn) with high threshold of similar size (ht) 
    allcomps = pd.DataFrame(compres)
    allcomps.to_csv(workdir+"supplementary_information/cn_ht_comparison.csv", sep=';', index=True, encoding='utf-8')
    
    allcomps = allcomps.T
    fig, axs = plt.subplots(2,1, figsize=(10,5), sharex=True, sharey=False)
    axs[0].plot(allcomps.index, allcomps["cn_cc"], linestyle="-", color="gray", marker=".")
    axs[0].plot(allcomps.index, allcomps["ht_cc"], linestyle="--", color="gray", marker=".")
    axs[0].plot(allcomps.index, allcomps["cn_pos_linkage"], linestyle="-", color="red", marker=".")
    axs[0].plot(allcomps.index, allcomps["ht_pos_linkage"], linestyle="--", color="red", marker=".")
    axs[0].plot(allcomps.index, allcomps["cn_neg_linkage"], linestyle="-", color="blue", marker=".")
    axs[0].plot(allcomps.index, allcomps["ht_neg_linkage"], linestyle="--", color="blue", marker=".")
    axs[1].plot(allcomps.index, allcomps["cn_modularity"], linestyle="-", color="orange", marker="o")
    axs[1].plot(allcomps.index, allcomps["ht_modularity"], linestyle="--", color="orange", marker="o")
    ax2 = axs[1].twinx()
    ax2.plot(allcomps.index, allcomps["cn_lcc"], linestyle="-", color="green", marker="o")
    ax2.plot(allcomps.index, allcomps["ht_lcc"], linestyle="--", color="green", marker="o")
    
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color="red", lw=4),
                    Line2D([0], [0], color="blue", lw=4),
                    Line2D([0], [0], color="gray", lw=4),
                    Line2D([0], [0], color="black", lw=1, linestyle="--"),
                    Line2D([0], [0], color="black", lw=1),]
    
    axs[0].legend(custom_lines, ['Interglacial Linkage', 'Glacial Linkage', 'Clustering Coefficient', 'Spiec-Easi Sparsification', 'Consensus Sparsification'])
    
    axs[0].set_ylabel("Clustering")
    axs[1].set_ylabel("Modularity")
    ax2.set_ylabel("LCC Size")
    axs[1].set_xlabel("Consensus Threshold")
    
    
    custom_lines = [Line2D([0], [0], color="green", lw=4),
                    Line2D([0], [0], color="orange", lw=4),
                    Line2D([0], [0], color="black", lw=1, linestyle="--"),
                    Line2D([0], [0], color="black", lw=1),]
    
    axs[1].legend(custom_lines, ['LCC Size', 'Modularity', 'Spiec-Easi Sparsification', 'Consensus Sparsification'], loc="center left")
    fig.savefig(workdir+'plots/compare_BN_HT.png', format='png', dpi=1200)

plt.show()

plotNet(tpos, tneg)
exportToGephi(thisG, df_families, tpos, tneg, thisMods, path+"gephi_spieceasi_05.gml", trophicLevels)

flowMetric = "A/DC"
proxvals = ngrip["temp"].values
ages = list(abunds.index)
sealevel = pd.read_csv("C:/Users/vdinkel/Desktop/Manuscript/input/Grant_age_RSL_col_EF.txt", delimiter="\t")
sealevel["(ka BP)"] = sealevel["(ka BP)"]* 1000
ageindeces = [np.abs(sealevel["(ka BP)"] - k).argmin() for k in ages] # find closest ages from sealevel df
sealevels = sealevel.iloc[ageindeces]["(m)"]


#grant_EF[("(ka BP)")] = grant_EF[("(ka BP)")]*1000
#newageindeces = [np.abs(grant_EF[("(ka BP)")] - k).argmin() for k in ages] # find closest ages from sealevel df
#sealevels = grant_EF.iloc[newageindeces]["(m)"].values

#excludeFams = ["Bovichtidae", "Hapalosiphonaceae", "Nostocaceae", "Salpingoecidae", "Microcoleaceae", "Oscillatoriaceae"]
#for exFam in excludeFams:
#    lcc.remove(exFam)

#for com_i in range(0,4):
flows, ascss, all_flows, all_mods_flows, lcc_lccs, lcc_maturity = flowIndeces(abunds, proxvals, trDict2, thisG.subgraph(lcc), flowMetric) #lcc #mods2["comms"][3]

def corrMaturity(maturity, proxvals, sealevels):
    print ("Correlation ---- TEMP")
    for key in maturity.keys():
        rt, pt = scipy.stats.spearmanr(proxvals, maturity[key])
        print (key, " = ", "r(", np.round(rt,2), "), p(" , np.round(pt,2),")")
    print ("Correlation ---- SEALEVEL")
    for key in maturity.keys():
        r, p = scipy.stats.spearmanr(sealevels, maturity[key])
        print (key, " = ", "r(", np.round(r,2), "), p(" , np.round(p,2),")")
    return {"temp": [np.round(rt,2) , np.round(pt,2)], "sealevel": [np.round(r,2), np.round(p,2)]}

def plotMaturity(ax1, maturity, title):
    # plots maturity for each module spererately
    ax2 = ax1.twinx()
    ax1.set_title(title)
    for key in maturity.keys():
        if key != "T" and key != "A/DC":
            ax1.plot(list(abunds.index), maturity[key], label=key, linestyle='dashed', linewidth=1,)
        else:
            ax2.plot(list(abunds.index), maturity[key], label=key)
            #ax2.set_ylabel("XY", color="black")
    
        #ax1.set_ylabel(flowMetric, color="black")
    ax1.legend()
    ax2.legend()
    
    ax1.get_xaxis().set_visible(False)
    ax2.get_xaxis().set_visible(False)
    
    ax1.spines['top'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)

def plotMaturity2(ax1, maturity, title):
    #plots maturity for each parameter (all modules together for each parameter)
    #ax1.set_title(title)
    labels = ["LCC", "Module 1","Module 3", "Module 6"]
        
    if title == "T":
        ytitle = "Total System Throughput\nTST"
    if title == "A/DC":
        ytitle = "Relative Ascendency\nA/DC"
    if title == "A":
        ytitle = "Ascendency\nA"
    if title == "DC":
        ytitle = "Development Capacity\nDC"
    if title == "H":
        ytitle = "Complexity\nH"
    
    i = 0
    sealevels_centered = (sealevels - np.mean(sealevels))/np.std(sealevels)

    correlations = []
    print("---------------------")
    print(title, "\t", "ngrip_R", "\t", "ngrip_P", "\t", "sealevel_R", "\t", "sealevel_P")
    for label in labels:
        ngrip_correlated = False
        sealevel_correlated = False
        linestyle='dashed'
        linewidth = 1
        corradd = ""
        # correlation with ngrip
        ngrip_corr = scipy.stats.spearmanr(ngrip["temp"].values, maturity[i][title])
        sealevels_corr = scipy.stats.spearmanr(sealevels_centered, maturity[i][title])
        if ngrip_corr[1] <= 0.05:
            linestyle = 'solid'
            linewidth = 2
            ngrip_correlated = True
            corradd += "p_ngrip*" 
            if ngrip_corr[1] <= 0.01:
                corradd += "*"
            if ngrip_corr[1] <= 0.001:
                corradd += "*"
        if sealevels_corr[1] <= 0.05:
            linestyle = 'solid'
            linewidth = 2
            sealevel_correlated = True
            if ngrip_correlated:
                corradd += ", "
            corradd += "p_rsl*" 
            if sealevels_corr[1] <= 0.01:
                corradd += "*"
            if sealevels_corr[1] <= 0.001:
                corradd += "*"
        #import pdb; pdb.set_trace()
        if corradd != "":
            corradd = " ("+corradd+")"
        
        print(label, "\t", round(ngrip_corr[0],2), "\t", round(ngrip_corr[1],2), "\t", round(sealevels_corr[0],2), "\t", round(sealevels_corr[1],2))
        correlations.append(
            {"metric": title,
             "module": label, 
             "ngrip_R": ngrip_corr[0],
             "ngrip_p": ngrip_corr[1],
             "sealevel_R": sealevels_corr[0],
             "sealevel_P": sealevels_corr[1]
             })
        # multiple linear regression maturity = ngrip + sealevel
        #import scipy.stats as ss
        #dd = {'y': maturity[i][title][1], 'ngrip': ss.rankdata(ngrip["temp"].values), 'sealevels': ss.rankdata(sealevels_centered)}
        #df_dd = pd.DataFrame(data=dd)
        #lm = smf.ols(formula='y ~ ngrip + sealevels', data=df_dd).fit()
        
        #print(title, " ", label)
        #print(lm.params)
        #print(lm.pvalues)
        
        #ax1.vlines(grid_ages, min(maturity[i][title]), max(maturity[i][title]), color="lightgray")
        ax1.plot(list(abunds.index), maturity[i][title], label=label+corradd, linestyle=linestyle, linewidth=linewidth, marker = "o", color = colmap[i])
        #if ngrip_correlated:
        #    text = ax1.text(list(abunds.index)[-1]-1000, maturity[i][title][len(maturity[i][title]) -1 ], "*", color="red", size=20)
        #if sealevel_correlated:
        #    text = ax1.text(list(abunds.index)[-1]-5000, maturity[i][title][len(maturity[i][title]) -1 ], "*", color="blue", size=20)
        
        ax1.set_ylabel(ytitle)
        i += 1
    
    ax1.legend(bbox_to_anchor=(1.0, .5))
    ax1.spines['top'].set_visible(False)
    if title != "DC":
        ax1.get_xaxis().set_visible(False)
        ax1.spines['bottom'].set_visible(False)
    else:
        ax1.set_xlabel("Years BP")
    return correlations

grid_ages= [0, 20000, 40000, 60000, 80000, 100000, 120000]
parameters = ["A/DC", "A", "DC"] #, "AMI", "H" "T",
fig, axs = plt.subplots(len(parameters)+1,1, figsize=(10,10), sharex=True, sharey=False)
#axs[0].vlines(grid_ages, min( ngrip["temp"].values), max( ngrip["temp"].values), color="lightgray")
axs[0].plot(ngrip["age"].values, ngrip["temp"].values, color = "black", label="NGRIP", alpha =0.8, linewidth = 2, marker="o")
axs[0].set_ylabel("Temperature\nNGRIP")
ax2 = axs[0].twinx()
#axs[1].vlines(grid_ages, min(sealevels), max(sealevels), color="lightgray")
ax2.plot(ngrip["age"].values, sealevels, label="RSL", color = "black", linestyle="--", linewidth = 2, alpha =0.8, marker="o")
ax2.set_ylabel("Relative Sea Level\nRSL")
axs[0].legend(loc="center", bbox_to_anchor=(.5,  .7))
ax2.legend(loc="center", bbox_to_anchor=(.5, 0.55))

axs[0].get_xaxis().set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[0].spines['bottom'].set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].spines['top'].set_visible(False)
axs[1].spines['bottom'].set_visible(False)
ax2.get_xaxis().set_visible(False)
ax2.spines['bottom'].set_visible(False)

mod_flows = []
for com_i in [0,2,5]:
    subgraph = thisG.subgraph(mods2["comms"][com_i])
    #if com_i == 6:
    #    subgraph = thisG.subgraph(list(mods2["comms"][2]) + list(mods2["comms"][4]) + list(mods2["comms"][6]))
    com_fams = list(mods2["comms"][com_i])
    flows, ascss, all_flows, all_mods_flows, lccs, mod_maturity = flowIndeces(abunds, proxvals, trDict2, subgraph, flowMetric) #lcc #mods2["comms"][3]
    mod_flows.append(mod_maturity)
    # supplements: module maturity parameters
    df_maturity = pd.DataFrame(mod_maturity)
    df_maturity.index = ages
    df_maturity.index.name = "years_BP"
    df_maturity.to_csv(workdir+"supplementary_information/maturity_module"+str(com_i+1)+".csv", sep=',', index=True, encoding='utf-8')

# supplements: lcc maturity parameters
df_maturity = pd.DataFrame(lcc_maturity)
df_maturity.index = ages
df_maturity.index.name = "years_BP"
df_maturity.to_csv(workdir+"supplementary_information/maturity_module"+str(com_i+1)+".csv", sep=',', index=True, encoding='utf-8')


colmap= ["#BBBBBB", "#CC3311", "#EE7733", "#0077BB"]
k = 1
corrs = []
for par in parameters:
    maturity_corrs = plotMaturity2(axs[k], [lcc_maturity, mod_flows[0], mod_flows[1], mod_flows[2]], par)
    corrs.append(maturity_corrs)
    k+=1

#scipy.stats.spearmanr(proxvals, mod_sum_flows[0])
fig.tight_layout()
fig.subplots_adjust(hspace=0.06)
fig.savefig(workdir+'plots/maturity.png', format='png', dpi=1200)
df_corrs = pd.DataFrame(corrs[0]).append(pd.DataFrame(corrs[1])).append(pd.DataFrame(corrs[2]))
df_corrs.to_csv(workdir+"supplementary_information/maturity_corrs.csv", sep=',', index=True, encoding='utf-8')
plt.show()

#import pdb; pdb.set_trace()

def getModuleProxies(module):
    thisMod = list(module)
    modAbunds = abunds[thisMod].sum(axis=1).values
    
    biogenic = biogenic_all["X..biogenic"].values
    silic = biogenic_all["X..siliciclastics"].values
    baal = xrf_all['ln.Ba.Al.'].values
    thisngrip = ngrip["temp"].values
    
    ret = {}
    r, p = scipy.stats.spearmanr(modAbunds, biogenic)
    if p <= 0.05:
        ret["biogenic"] = int(r*100)/100
    r, p = scipy.stats.spearmanr(modAbunds, silic)
    if p <= 0.05:
        ret["silic"] = int(r*100)/100
    r, p = scipy.stats.spearmanr(modAbunds, baal)
    if p <= 0.05:
        ret["ba/al"] = int(r*100)/100
    r, p = scipy.stats.spearmanr(modAbunds, thisngrip)
    if p <= 0.05:
        ret["ngrip"] = int(r*100)/100
        
    ret["posneg"] = {"biogenic": {"pos": [], "neg": []}, "silic": {"pos": [], "neg": []}, "ba/al": {"pos": [], "neg": []}, "ngrip": {"pos": [], "neg": []}}
    for tax in thisMod:
        taxvals = abunds[thisMod[0]].values
        # ngrip
        r, p = scipy.stats.pearsonr((taxvals-np.mean(taxvals))/np.std(taxvals), thisngrip)
        if p < 0.9 and abs(r)>0.2:
            if r > 0:
                ret["posneg"]["ngrip"]["pos"].append(tax)
            else:
                ret["posneg"]["ngrip"]["neg"].append(tax)
        # biogenic
        r, p = scipy.stats.pearsonr((taxvals-np.mean(taxvals))/np.std(taxvals), biogenic)
        if p < 0.9 and abs(r)>0.2:
            if r > 0:
                ret["posneg"]["biogenic"]["pos"].append(tax)
            else:
                ret["posneg"]["biogenic"]["neg"].append(tax)
        # silic
        r, p = scipy.stats.pearsonr((taxvals-np.mean(taxvals))/np.std(taxvals), silic)
        if p < 0.9 and abs(r)>0.2:
            if r > 0:
                ret["posneg"]["silic"]["pos"].append(tax)
            else:
                ret["posneg"]["silic"]["neg"].append(tax)
        # baal
        r, p = scipy.stats.pearsonr((taxvals-np.mean(taxvals))/np.std(taxvals), baal)
        if p < 0.9 and abs(r)>0.2:
            if r > 0:
                ret["posneg"]["ba/al"]["pos"].append(tax)
            else:
                ret["posneg"]["ba/al"]["neg"].append(tax)
        #import pdb; pdb.set_trace()
    return ret

def plotModule(comm, ps, ns):
    
    modProxies = getModuleProxies(module)
    #ps = modProxies["posneg"]["biogenic"]["pos"]
    #ns = modProxies["posneg"]["biogenic"]["neg"]
    
    x_offset = 300
    y_offset = 1000
    comms = [list(k) if len(list(k))>1 else [] for k in thisMods['comms']]
    
    plt.figure(figsize=[18,10], dpi = 300)
    pos = nx.spring_layout(thisG.subgraph(comm), center = [0,0], seed = 6, scale=500) #4
    #nx.draw_networkx_nodes(G_masked, pos)
    #plt.show()

    trColors = ["#1D8348", "#D4AC0D", "#A04000", "#884EA0"]
    trHeights = 500
    
    i = 0
    flipFlop = 1
    colorBy = "ngrip"#"ngrip"
    
    plt.axhspan(-500, 500, facecolor='0.2', alpha=0.1)
    plt.axhspan(1500, 2500, facecolor='0.2', alpha=0.1)

    maxNodeSize = 300
    commPadding = 100
    prevStartX = 0
    prevEndX = 0
    allTRs = {"1": [], "2": [], "3": [], "4": []}
    
    # comm
    colors = []
    randCol = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])][0]
    trTaxa = {}
    trLevels = {"1": [], "2": [], "3": [], "4": []}
    trGrps = {"1": [], "2": [], "3": [], "4": []}
    
    # get trophic levels of community
    for k in comm:
        trGrp = df_families[df_families[0] == k][1].values[0]
        trLevel = trophicLevels[trophicLevels[0] == trGrp][1].values[0]
        trTaxa[k] = trLevel
        trLevels[str(trLevel)].append(k)
        trGrps[str(trLevel)].append(trGrp)
        allTRs[str(trLevel)].append(k)
    trCounts = Counter([trTaxa[k] for k in trTaxa.keys()])
    maxCounts = trCounts.most_common()[0][1]
    
    node_sizes = dict(np.sum(abunds[comm])*maxNodeSize)
    for key in node_sizes.keys():
        node_sizes[key] = max(min(node_sizes[key], maxNodeSize), maxNodeSize/6)
    moduleSizeSum = sum([node_sizes[k] for k in trLevels[str(trCounts.most_common()[0][0])]])
    
    commStartX = 0
    commEndX = 1000
    commWidth = 1000
    prevStartX = 0
    prevEndX = 1000
    
    for posKey in pos.keys():
        pos[posKey][0] = pos[posKey][0] + 500
    
    for tr in trLevels.keys():
        y_base =  y_offset*(int(tr)-1)
        y_upper = y_base + trHeights - 100
        y_lower = y_base - trHeights + 100
        
        taxaDict = {}
        for tax in trLevels[tr]:
            taxaDict[tax] = pos[tax][1]
        taxaDict = {k: v for k, v in sorted(taxaDict.items(), key=lambda item: item[1])} # sorted
        
        try:
            spaceBetween = abs(y_upper - y_lower) / (len(trLevels[tr])-1)
        except:
            spaceBetween = trHeights
        for tax in taxaDict.keys():
            pos[tax][1] = int(y_lower+(spaceBetween * list(taxaDict.keys()).index(tax)))
        
        #import pdb; pdb.set_trace()
        
    for k in comm:
        
        trLevel = trTaxa[k]
        if colorBy == "module":
            colors.append(randCol)  
        if colorBy == "trophic":
            colors.append(trColors[trLevel-1])  
        if colorBy == "ngrip":
            if k in ps:
                colors.append(np.array([0.9, 0.0, 0.0, 0.85]))  
            elif k in ns:
                colors.append(np.array([0.0, 0.0, 0.7, 0.85]))
            else:
                colors.append(np.array([0.7, 0.7, 0.7, 0.7]))
        if colorBy == "bio":
            if k in ps:
                colors.append(np.array([0.0, 0.39, 0.0, 0.7]))  
            elif k in ns:
                colors.append(np.array([0.58, 0.0, 0.83, 0.7]))
            else:
                colors.append(np.array([0.7, 0.7, 0.7, 0.7]))

    nx.draw_networkx_nodes(thisG, pos, nodelist=comm, node_color=colors, node_size=list(node_sizes.values()))
    
    for node in comm:
        labelSize = max((node_sizes[node] / maxNodeSize)*13, 7)
        plt.text(pos[node][0], pos[node][1], node, fontsize=labelSize, ha='center', va='center')
    i+=1
    flipFlop *= -1
    
    edgecols = []
    alphas = []
    for edge in thisG.subgraph(comm).edges:
        baseCol = np.array([0.7, 0.7, 0.7, 0.7])
        if colorBy == "ngrip":
            if edge[0] in ps or edge[1] in ps:
                baseCol = np.array([0.9, 0.0, 0.0, 0.7])
            if edge[0] in ns or edge[1] in ns:
                baseCol = np.array([0.0, 0.0, 0.9, 0.7])
            if (edge[0] in ps and edge[1] in ns) or (edge[0] in ns and edge[1] in ps):
                baseCol = np.array([60/255, 0.0, 100/255, 0.7])
        edgecols.append(baseCol)
        alphas.append(multiplex_A[edge[0]][edge[1]])
    
    nx.draw_networkx_edges(thisG.subgraph(comm), pos, edge_color = edgecols, width=np.array(alphas)*2, alpha=np.array(alphas))
    plt.xlim(0-150, commEndX+150)
    plt.ylim(-500, 3500)
    import time
    obj = time.gmtime(0)
    curr_time = round(time.time()*1000)
    plt.savefig("C:/Users/vdinkel/Desktop/PhD/Presentations/networkVisialization/"+str(curr_time)+".png")
    
    assocCover = getModuleAssociationCoverage(module)
    modTemp = getModuleTemp(module)
    
    
    txtStartY = 3300
    txtYOffset = 150
    
    plt.text(0, txtStartY, "Size: "+str(len(module)), bbox=dict(facecolor='white', alpha=0.5))
    plt.text(0, txtStartY-txtYOffset, "MTL: "+str(getModuleMTL(module)[0]), bbox=dict(facecolor='white', alpha=0.5))
    plt.text(0, txtStartY-txtYOffset*2, "MTL (layered): [1]: "+str(getModuleMTL(module)[1][0])+" | [2]: "+str(getModuleMTL(module)[1][1])+" | [3]: "+str(getModuleMTL(module)[1][2])+" | [4]: "+str(getModuleMTL(module)[1][3]), bbox=dict(facecolor='white', alpha=0.5))
    plt.text(0, txtStartY-txtYOffset*3, "Association Coverage [Direct]: "+str(assocCover[1]["direct"]*100)+"% | [Indirect] "+str(assocCover[1]["indirect"]*100)+"%", bbox=dict(facecolor='white', alpha=0.5))
    plt.text(0, txtStartY-txtYOffset*4, "Ngrip [Interglacial]: "+str(modTemp["pos"])+" | [Glacial] "+str(modTemp["neg"]), bbox=dict(facecolor='white', alpha=0.5))
    
    i = 5
    if 'ba/al' in modProxies.keys():
        plt.text(0, txtStartY-txtYOffset*i, "Ba/Al: "+str(modProxies["ba/al"]), bbox=dict(facecolor='white', alpha=0.5))
        i += 1
    if 'biogenic' in modProxies.keys():
        plt.text(0, txtStartY-txtYOffset*i, "Biogenic: "+str(modProxies["biogenic"]), bbox=dict(facecolor='white', alpha=0.5))
        i += 1
    if 'silic' in modProxies.keys():
        plt.text(0, txtStartY-txtYOffset*i, "Silic: "+str(modProxies["silic"]), bbox=dict(facecolor='white', alpha=0.5))
        i += 1
    if 'ngrip' in modProxies.keys():
        plt.text(0, txtStartY-txtYOffset*i, "Ngrip: "+str(modProxies["ngrip"]), bbox=dict(facecolor='white', alpha=0.5))
        i += 1
    
    plt.show()
    plt.close()