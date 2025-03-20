# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 15:31:44 2024

@author: vdinkel
"""

import networkx as nx
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import pygraphviz
from adjustText import adjust_text

def calcAsc(G):

    T =sum([k[2]["weight"] for k in G.edges(data=True)]) # sum of all edge weights
    AMI = 0.0
    Hc = 0.0
    H = 0.0
    DC = 0.0
    
    edgevals = {}
    
    for edge in G.edges(data=True):
        Tj = 0.0
        Ti = 0.0
        for inedge in G.in_edges(edge[1], data=True):
            Tj += inedge[2]["weight"]
        for outnedge in G.out_edges(edge[0], data=True):
            Ti += outnedge[2]["weight"]   
        Tij = edge[2]["weight"]
        
        thisAMI = Tij/T * np.log2((Tij * T)/(Ti * Tj))
        thisHc = -1 * (Tij/T * np.log2((Tij*Tij)/T))
        thisH = thisAMI + thisHc
        thisDC = T * thisH
        
        AMI += thisAMI
        Hc += thisHc
        H += thisH
        DC += thisDC
        
        edgevals[(edge[0],edge[1])] = {"AMI": thisAMI, "Hc": thisHc, "H": thisH, "DC": thisDC, "Tij": Tij, "T": T}

    return T, AMI, Hc, H, DC, edgevals

def getSimG(G, vabunds):
    diG = nx.DiGraph()
    nodes = list(G.nodes())
    for node in nodes:
        if vabunds[node]>0.0:
            diG.add_node(node)
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
        diG[edge[0]][edge[1]]["weight"] = vabunds[edge[0]] / degs[edge[0]]
    return diG

def flowIndeces(flowG, abunds):
    # proxievals = ngrip["temp"].values
    AMIs = []
    Ts = []
    Hs = []

    for k in range(0,len(abunds)):
        vabunds = abunds.iloc[[k]]
        simG = getSimG(flowG, vabunds)
        
        T, AMI, Hc, H, DC =  calcAsc(simG)
        AMIs.append(AMI)
        Ts.append(T)
        Hs.append(H)
    
    As = np.array(AMIs) * np.array(Ts)
    DCs = np.array(Ts) * np.array(Hs)
    '''
    print("-----------")
    print("AMI:  -> ", scipy.stats.pearsonr(AMIs, proxievals))
    print("T:  -> ", scipy.stats.pearsonr(Ts, proxievals))
    print("A: -> ", scipy.stats.pearsonr(As, proxievals))
    print("H: -> ", scipy.stats.pearsonr(Hs, proxievals))
    print("DC:  -> ", scipy.stats.pearsonr(DCs, proxievals))
    print("A/DC:  -> ", scipy.stats.pearsonr(As/DCs, proxievals))
    '''
    
    maturity = {"AMI": AMIs, "T": Ts, "A": As, "H": Hs, "DC" : DCs, "A/DC": As/DCs }
    
    return maturity

def cluster_module(df_module, ages):
    # Scale the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df_module)
    scaled_df = pd.DataFrame(scaled_data, columns=df_module.columns)
    
    # Apply K-Means
    kmeans = KMeans(n_clusters=2, random_state=12)  # Set random_state for reproducibility
    kmeans.fit(scaled_data)
    
    # Add cluster labels to the original DataFrame
    df_module['Cluster'] = kmeans.labels_
    df_module['years_BP'] = ages
    return df_module

def exportToGephi(matret, fileName, troph, c):
    
    G = matret["simG"+c]
    abunds =  matret["abunds_c"+c]
    edgevals = matret["edgevals_c"+c]
    
    f = open(fileName, "w+")
    f.write("graph\n[\n")
    # NODES
    
    degrees = dict(G.degree)
    nodes = list(set([k if degrees[k] > 0 else None for k in degrees.keys()]))
    for node in nodes:
        if node != None:
            # !!! MODULE !!!
            #moduleID = getNodeModuleID(node, modules)
            f.write("node\n[\n")
            f.write('id "'+node+'"\n')
            f.write('label "'+node+'"\n')
            f.write('trophic_level '+str(troph[troph["Family"] == node]["Level"].values[0])+'\n')
            f.write("abund "+str(abunds[node])+'\n')
            f.write("]\n")
    
    # EDGES
    for edge in G.edges:
       
        f.write("edge\n[\n")
        f.write('source "'+edge[0]+'"\n')
        f.write('target "'+edge[1]+'"\n')

        weight = G.get_edge_data(edge[0], edge[1])["weight"]
        f.write('label "'+str(int(weight))+'"\n')
        
        edge_A = edgevals[edge]["AMI"] * edgevals[edge]["T"]
        edge_ADC = edge_A / edgevals[edge]["DC"]
        
        
        f.write('weight '+str(weight)+'\n')
        f.write('AMI '+str(edgevals[edge]["AMI"])+'\n')
        f.write('Hc '+str(edgevals[edge]["Hc"])+'\n')
        f.write('H '+str(edgevals[edge]["H"])+'\n')
        f.write('DC '+str(edgevals[edge]["DC"])+'\n')
        f.write('Tij '+str(edgevals[edge]["Tij"])+'\n')
        try:
            	f.write('diff '+str(edgevals[edge]["diff"])+'\n')
        except:
            f.write('diff '+str(1.0)+'\n')
        f.write('A '+str(edge_A)+'\n')
        f.write('A/DC '+str(edge_ADC)+'\n')

        f.write("]\n")
    
    f.write("]")
    f.close()
    
def z_score(values, x):
    values = np.array(values)
    z = (x - np.mean(values)) / np.std(values)
    return round(z,2)

# LOAD TROPHIC META
# ADD TROPHIC
netsdir = "C:/Users/vdinkel/Desktop/Manuscript2/scripts/man1/"
path = "C:/Users/vdinkel/Desktop/Manuscript/submission_Nature_Communications/revision1/script/"
trophic = pd.read_csv(netsdir+"trophic_levels.csv", sep=",", header=None, index_col=False, names=['Group', 'Level'])
fam_trophic = pd.read_csv(netsdir+"families.csv", sep=";", header=None, index_col=False, names=['Family', 'Group'])
# Merge Globo & Troph
troph = pd.merge(fam_trophic, trophic, left_on="Group", right_on="Group", how='left')
# ////ADD TROPHIC
ages = pd.read_csv(path+"Table_S2_maturity_moduleIG1.csv")["years_BP"]

# NGRIP
ngrip = pd.read_csv("C:/Users/vdinkel/Desktop/Manuscript/submission_Nature_Communications/code/consensus_network/output/KL-77_ngrip.csv", sep=";")

abunds = pd.read_table("C:/Users/vdinkel/Desktop/Manuscript/output/F_KL-77_rel.csv", sep=";", header=0)
G = nx.read_gml(path+"cn_spieceasi_05.gml")

nodes_with_module_ig1 = [n for n, attr in G.nodes(data=True) if attr.get("module") == "0"]
nodes_with_module_ig2 = [n for n, attr in G.nodes(data=True) if attr.get("module") == "5"]
nodes_with_module_glc =[n for n, attr in G.nodes(data=True) if attr.get("module") == "2"]

G_ig1 = G.subgraph(nodes_with_module_ig1)
G_ig2 = G.subgraph(nodes_with_module_ig2)
G_glc = G.subgraph(nodes_with_module_glc)

maturity_ig1 = pd.read_csv(path+"Table_S2_maturity_moduleIG1.csv")[["A/DC", "A", "DC"]]
maturity_ig2 = pd.read_csv(path+"Table_S4_maturity_moduleIG2.csv")[["A/DC"]]
maturity_glc = pd.read_csv(path+"Table_S3_maturity_moduleGLC.csv")[["A/DC"]]

maturity_ig1["ngrip"] = (ngrip["temp"].values).reshape(-1, 1)
maturity_ig2["ngrip"] = (ngrip["temp"].values).reshape(-1, 1)
maturity_glc["ngrip"] = (ngrip["temp"].values).reshape(-1, 1)

maturity_ig1 = cluster_module(maturity_ig1, ages)
maturity_ig2 = cluster_module(maturity_ig2, ages)
maturity_glc = cluster_module(maturity_glc, ages)

def getEdgeMaturities(G, nodes, maturity, abunds):
    c0 = maturity[maturity["Cluster"] == 0]["years_BP"].values
    c0_abunds = np.mean(abunds.loc[c0][nodes])
    
    simG0 = getSimG(G, c0_abunds)
    c0_T, c0_AMI, Hc, H, c0_DC, edgevals_c0 =  calcAsc(simG0)
    
    c1 = maturity[maturity["Cluster"] == 1]["years_BP"].values
    c1_abunds = np.mean(abunds.loc[c1][nodes])
    
    simG1 = getSimG(G, c1_abunds)
    c1_T, c1_AMI, Hc, H, c1_DC, edgevals_c1 =  calcAsc(simG1)
    
    #c0_sum_DC = sum([edgevals_c0[k]["DC"] for k in edgevals_c0.keys()])
    #c1_sum_DC = sum([edgevals_c1[k]["DC"] for k in edgevals_c1.keys()])
    
    diffkey = "Tij"
    for edge in edgevals_c0.keys():
        if edge not in edgevals_c1.keys():
            edgevals_c0[edge]["diff"] = 10
            edgevals_c0[edge]["c0_val"] = 0
            edgevals_c0[edge]["c1_val"] = 0
        else:
            diffval_c0 = (edgevals_c0[edge]["AMI"] * edgevals_c0[edge]["T"]) / c0_DC #edgevals_c0[edge]["DC"] #edgevals_c0[edge][diffkey]#
            diffval_c1 = (edgevals_c1[edge]["AMI"] * edgevals_c1[edge]["T"]) / c1_DC #edgevals_c1[edge]["DC"] #edgevals_c1[edge][diffkey]#
            edgevals_c0[edge]["diff"] = diffval_c0 - diffval_c1
            edgevals_c0[edge]["c0_perc"] = (diffval_c0 / ((c0_AMI * c0_T)/c0_DC))
            edgevals_c0[edge]["c1_perc"] = (diffval_c1 / ((c1_AMI * c1_T)/c1_DC))
            edgevals_c0[edge]["c0_val"] = diffval_c0
            edgevals_c0[edge]["c1_val"] = diffval_c1
    
    for edge in edgevals_c1.keys():
        if edge not in edgevals_c0.keys():
            edgevals_c1[edge]["diff"] = 10
            edgevals_c1[edge]["c0_val"] = 0
            edgevals_c1[edge]["c1_val"] = 0
        else:
            diffval_c0 = (edgevals_c0[edge]["AMI"] * edgevals_c0[edge]["T"]) /c0_DC # edgevals_c0[edge]["DC"] #edgevals_c0[edge][diffkey]#
            diffval_c1 = (edgevals_c1[edge]["AMI"] * edgevals_c1[edge]["T"]) / c1_DC # edgevals_c1[edge]["DC"] #edgevals_c1[edge][diffkey]#
            edgevals_c1[edge]["diff"] = diffval_c1 - diffval_c0
            edgevals_c1[edge]["c0_perc"] = (diffval_c0 / ((c0_AMI * c0_T)/c0_DC))
            edgevals_c1[edge]["c1_perc"] = (diffval_c1 / ((c1_AMI * c1_T)/c1_DC))
            edgevals_c1[edge]["c0_val"] = diffval_c0
            edgevals_c1[edge]["c1_val"] = diffval_c1

    return {"simG0": simG0, "edgevals_c0": edgevals_c0, "abunds_c0": c0_abunds, "c0_ADC": ((c0_AMI*c0_T)/c0_DC), "simG1": simG1, "edgevals_c1": edgevals_c1, "abunds_c1": c1_abunds, "c1_ADC": ((c1_AMI*c1_T)/c1_DC)}


ig1_matret = getEdgeMaturities(G_ig1, nodes_with_module_ig1, maturity_ig1, abunds)
exportToGephi(ig1_matret, path+"/graphs/simIG1_c0.gml", troph, c="0")
exportToGephi(ig1_matret, path+"/graphs/simIG1_c1.gml", troph, c="1")

ig2_matret = getEdgeMaturities(G_ig2, nodes_with_module_ig2, maturity_ig2, abunds)
exportToGephi(ig2_matret, path+"/graphs/simIG2_c0.gml", troph, c="0")
exportToGephi(ig2_matret, path+"/graphs/simIG2_c1.gml", troph, c="1")

glc_matret = getEdgeMaturities(G_glc, nodes_with_module_glc, maturity_glc, abunds)
exportToGephi(glc_matret, path+"/graphs/simGLC_c0.gml", troph, c="0")
exportToGephi(glc_matret, path+"/graphs/simGLC_c1.gml", troph, c="1")

# PLOT NGRIP CLUSTERS
fig, axs = plt.subplots(4, 1, sharex=True, tight_layout=True, figsize=(10, 4))
plt.xticks(rotation=75)

ig1_counts = maturity_ig1["Cluster"].values
ig2_counts = maturity_ig2["Cluster"].values
glc_counts = maturity_glc["Cluster"].values

ngrip_col_warm = "#333333"
ngrip_col_cold = "#CCCCCC"
ngrip_colors = [ngrip_col_cold if k < 0 else ngrip_col_warm for k in ngrip["temp"].values]
ages_x = [str(k) for k in ages.values]

cluster_col_inc = "#333333"
cluster_col_dec = "#CCCCCC"
axs[0].bar(ages_x, ngrip["temp"].values, color=ngrip_colors, alpha=0.80) #, label=bar_labels, color=bar_colors
axs[1].bar(ages_x, [1 for k in ig1_counts], color=[cluster_col_inc if k == 1 else cluster_col_dec for k in ig1_counts]) # 0 = low, 1 = high
axs[2].bar(ages_x, [1 for k in ig2_counts], color=[cluster_col_inc if k == 1 else cluster_col_dec for k in ig2_counts]) # 0 = low, 1 = high
axs[3].bar(ages_x, [1 for k in glc_counts], color=[cluster_col_inc if k == 1 else cluster_col_dec for k in glc_counts]) # 0 = high, 1 = low

axs[0].set_ylabel("Standardized\nNGRIP δ¹⁸O")
axs[1].set_ylabel("IG1\nCluster")
axs[1].yaxis.set_ticks([])
axs[2].set_ylabel("IG2\nCluster")
axs[2].yaxis.set_ticks([])
axs[2].set_xlabel("kyrs BP")
plt.savefig("C:/Users/vdinkel/Desktop/Manuscript/submission_Nature_Communications/revision1/script/new_plots/S7_ngrip.png")
plt.show()

def plotMatret(matret, c_key, title): #c_key = "1" or "0", 
    node_colors = []
    trColors = ["#1D8348", "#D4AC0D", "#A04000", "#884EA0"]
    
    comms = {0: [], 1: [], 2: [], 3: []}
    for node in list(matret["simG"+c_key].nodes):
        trLevel = troph[troph["Family"] == node]["Level"].values[0]
        node_colors.append(trColors[trLevel-1])
        comms[trLevel-1].append(node)
    communities = [frozenset(comms[i]) for i in range(0,4)]
    
    plt.figure(figsize=(12,12))
    plt.title(title + " | Cluster = "+c_key)
    
    ''' NEW
    '''
    G = matret["simG"+c_key]
    supergraph = nx.cycle_graph(4)
    superpos = nx.spring_layout(supergraph, scale=2.4, seed=2)
    
    centers = list(superpos.values())
    centers[1][1] = -1
    centers[3][1] = -.8
    pos = {}
    ks = [0.5, 0.1, 3, 0.1]
    scales = [0.5, 0.5, 2.4, 0.5]
    i = 0
    for center, comm in zip(centers, communities):
        pos.update(nx.spring_layout(nx.subgraph(G, comm), center=center, scale = scales[i], k= ks[i], seed=1430))
        i+=1
    
    # Nodes colored by cluster
    if title == "IG1":
        maturity = maturity_ig1
    if title == "IG2":
        maturity = maturity_ig2
    if title == "GLC":
        maturity = maturity_glc
    kyrs = maturity[maturity["Cluster"] == int(c_key)]["years_BP"].values

    for nodes, clr in zip(communities, ("#1D8348", "#D4AC0D", "#A04000", "#884EA0")): #("#1D8348", "#D4AC0D", "#A04000", "#884EA0")
        #abund_sum = (abunds[node].loc[kyrs].values.sum()*1000000)
        nx.draw_networkx_nodes(G, pos=pos, nodelist=nodes, node_color=clr, node_size=700)
        #import pdb; pdb.set_trace()
    
    nx_labels = nx.draw_networkx_labels(G, pos, font_size=20, font_color="black")
    
    inc_col = "#333333" # orange ##D4A017
    dec_col = "#CCCCCC" # #1F77B4
    
    # TEST: GLACIAL / INTERGLACIAL : only bottom up / top down edges
    tmp_weights = [matret["edgevals_c"+c_key][edge]["diff"]*100 for edge in G.edges]
    tmp_colors = [inc_col if weight > 0 else dec_col for weight in tmp_weights]
    #tmp_indices = [i for i, color in enumerate(tmp_colors) if dec_col in color] # colors which are removed, so if dec_col then decreases are removed, only increases remain
    #edge_list = list(G.edges)
    #edges_to_remove = [edge_list[i] for i in tmp_indices]
    #G.remove_edges_from(edges_to_remove)      
        
    edges_to_remove = []
    for edge in G.edges:
        tr_src = troph[troph["Family"] == edge[0]]["Level"].values[0]
        tr_trg = troph[troph["Family"] == edge[1]]["Level"].values[0]
        if tr_src != tr_trg: # 
            edges_to_remove.append(edge)
    for edge in edges_to_remove:
        G.remove_edge(edge[0], edge[1])
    
    print("EDGES: ", len(G.edges))
    print("SUM: ", round(sum([matret["edgevals_c"+c_key][edge]["c"+c_key+"_val"] for edge in G.edges]),4))

    #print("C"+c_key+"_SUM: ", round(sum([matret["edgevals_c"+c_key][edge]["diff"] for edge in G.edges]),4))
    
    print("PERC_"+c_key+": ", round(sum([matret["edgevals_c"+c_key][edge]["c"+c_key+"_perc"] for edge in G.edges]),4))
    #print("PERC_1: ", round(sum([matret["edgevals_c1"][edge]["c1_perc"] for edge in G.edges]),4))
    
    print("C"+c_key+"_ADC: ", matret["c"+c_key+"_ADC"])
    #print("C1_ADC: ", matret["c1_ADC"])
    
    import pdb; pdb.set_trace()
    
    weights = [matret["edgevals_c"+c_key][edge]["diff"]*10000 for edge in G.edges]
    colors = [inc_col if weight > 0 else dec_col for weight in weights]
    weights = [max(.8, min(abs(weight), 8)) for weight in weights]
    #normalized_alphas = (np.array(weights) - min(weights)) / (max(weights) - min(weights))
    
    nx.draw_networkx_edges(G, pos, edge_color = colors, width=weights, arrowsize=30, connectionstyle="arc3,rad=0.2")

    plt.tight_layout()
    
    #adjust_text([nx_labels[key] for key in nx_labels.keys()])
    #plt.savefig("C:/Users/vdinkel/Desktop/Manuscript/submission_Nature_Communications/revision1/script/new_plots/S7_ig1_intraguild.png")
    plt.show()
    
    ''' BACKUP
    pos = nx.nx_agraph.graphviz_layout(matret["simG"+c_key], prog="neato", args='-Gnodesep=2')

    nx.draw_networkx_nodes(matret["simG"+c_key], pos, node_size=400, node_color=node_colors)
    nx.draw_networkx_labels(matret["simG"+c_key], pos, font_size=12, font_color="black")
    
    weights = [matret["edgevals_c"+c_key][edge]["diff"]*100 for edge in matret["simG"+c_key].edges]
    colors = ["green" if weight > 0 else "red" for weight in weights]
    weights = [max(.5, min(abs(weight), 10)) for weight in weights]
    
    nx.draw_networkx_edges(matret["simG"+c_key], pos, edge_color = colors, width=weights)
    '''
    plt.show()
 
#plotMatret(ig1_matret, "0", title = "IG1")
plotMatret(ig1_matret, "1", title = "IG1")

#plotMatret(ig2_matret, "0", title = "IG2")
#plotMatret(ig2_matret, "1" , title = "IG2")

#plotMatret(glc_matret, "0", title = "GLC")
#plotMatret(glc_matret, "1" , title = "GLC")