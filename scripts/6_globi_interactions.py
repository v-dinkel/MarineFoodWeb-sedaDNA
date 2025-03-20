# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 12:33:06 2024

@author: vdinkel
"""

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import json
import requests
import time
from collections import Counter

def read_config(filename):
    #read the config file of the project
    config = {}
    with open(filename, 'r') as file:
        for line in file:
            if line.strip() and not line.startswith("#"):  # Ignore empty lines or comments
                key, value = line.strip().split(":", 1)
                config[key.strip()] = value.strip()
    return config

def msum(A):
    return int(np.sum(np.sum(A)))

def AND(A, B):
    return B.where(A == 1, 0)

def loadConsensusEdgelist(path):
    df = pd.read_csv(path, sep=";")
    G = nx.from_pandas_edgelist(df, source='Source', target='Target')
    adj = nx.to_pandas_adjacency(G)
    return adj

def load_spiec_easi(direct):
    SE_A = pd.read_csv(direct+".csv", delimiter=";", header=0, index_col=0)
    SE_edgelist = pd.read_csv(direct+"_edgelist.csv", delimiter=";", header=0, index_col=0)
    
    SE_posneg = SE_A.copy()
    for idx in SE_edgelist.index:
        v1 = SE_edgelist.loc[idx]["v1"]
        v2 = SE_edgelist.loc[idx]["v2"]
        asso = SE_edgelist.loc[idx]["asso"]
        mult = 1
        if asso < 0:
            mult = -1
            SE_posneg.at[v1,v2] = SE_posneg.loc[v1][v2] * mult
    return SE_posneg

def getInteractionDB(net, globi):
    # this function returns for a given network (net) the reduced (occuring interaction) globi interactions database
    G_net = nx.from_pandas_adjacency(net)
    net_interactions = []
    i = 0
    for edge in G_net.edges:
        ret = globi[((globi["fam1"] == edge[0]) & (globi["fam2"] == edge[1])) | ((globi["fam1"] == edge[1]) & (globi["fam2"] == edge[0]))]
        if i == 0:
            net_interactions = ret.copy()
        if len(ret)>0:
            net_interactions = pd.concat([net_interactions, ret], ignore_index=True)
        i += 1
    net_interactions = net_interactions.drop_duplicates()
    return net_interactions

def mergeInteractions(G1, G2):
    # this function merges the interactions of an input graph G1 with GloBI graph G2
    # the weight of the interaction is the sum of occurences in G2
    
    # Normalize edges lexicographically so each edge is stored as (min, max)
    edges_G1 = {tuple(sorted(edge)) for edge in G1.edges()}
    edges_G2 = {tuple(sorted(edge)) for edge in G2.edges()}
    
    common_edges = edges_G1 & edges_G2
    
    weights = {edge: G2[edge[0]][edge[1]]["weight"] for edge in common_edges}
    weights = list(weights.values())
    levels = {edge: G2[edge[0]][edge[1]]["troph"] for edge in common_edges}
    levels = list(levels.values())
    
    troph_interacitons = {"1": 0, "2": 0, "3": 0, "4": 0}
    for i in range(0, len(weights)):
        troph_interacitons[str(levels[i])] += weights[i]
    
    for edge in G1.edges:
        G1.edges[edge[0], edge[1]]['weight'] = 0
        G1.edges[edge[0], edge[1]]['color'] = 'gray'
        G1.edges[edge[0], edge[1]]['troph'] = 0
        
    for edge in list(common_edges):
        G1.edges[edge[0], edge[1]]['weight'] = G2.edges[edge[0], edge[1]]["weight"]
        G1.edges[edge[0], edge[1]]['color'] = 'purple'
        G1.edges[edge[0], edge[1]]['troph'] = G2.edges[edge[0], edge[1]]["troph"]
        
    return G1, len(common_edges), troph_interacitons

def exportToGephi(G, fileName, troph, globi_list, cn_indirs, G_man1):
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
            f.write('module "'+G_man1.nodes[node]["module"]+'"\n')
            f.write('trophic_level "'+str(troph[troph["Family"] == node]["Level"].values[0])+'"\n')
            f.write("]\n")
    
    # EDGES
    for edge in G.edges:
        f.write("edge\n[\n")
        f.write('source "'+edge[0]+'"\n')
        f.write('target "'+edge[1]+'"\n')

        weight = G.get_edge_data(edge[0], edge[1])["weight"]
        if weight == 0:
            weight = 1
            f.write('label ""\n')
        else:
            f.write('label "'+str(int(weight))+'"\n')
            
        indir = 0
        indir_weight = 0
        if edge in cn_indirs.keys():
            indir = 1
            indir_weight = str(cn_indirs[edge])
            G.get_edge_data(edge[0], edge[1])["color"] = "orange"
        
        f.write('indir '+str(indir)+'\n')
        f.write('indir_weight '+str(indir_weight)+'\n')
            
        f.write('weight '+str(weight)+'\n')
        f.write('color '+G.get_edge_data(edge[0], edge[1])["color"]+'\n')

        f.write("]\n")
    
    f.write("]")
    f.close()
    
workdir = read_config("../config.txt")["workdir"]
outdir = workdir+"output/"
suppdir = workdir+"supplementary_information/"
indir = workdir+"input/"
plotdir = workdir+"plots/"
G_man1 = nx.read_gml(workdir+"output/cn_spieceasi_05.gml")

path = "C:/Users/vdinkel/Desktop/Manuscript/submission_Nature_Communications/revision1/script/globiDB_downloads/"
# Load GloBI file
load_df = pd.read_csv(indir+"db_allresult_final.csv")
predeats_df = load_df[load_df["interaction_type"].isin(["preysOn", "eats"])]
predeats_ = predeats_df[["fam1", "fam2", "source_taxon_name", "interaction_type", "target_taxon_name", "study_title"]]

# ADD TROPHIC
trophic = pd.read_csv(indir+"trophic_levels.csv", sep=",", header=None, index_col=False, names=['Group', 'Level'])
fam_trophic = pd.read_csv(indir+"families.csv", sep=";", header=None, index_col=False, names=['Family', 'Group'])

# Merge Globi & Troph
troph = pd.merge(fam_trophic, trophic, left_on="Group", right_on="Group", how='left')
predeats_ = pd.merge(predeats_, troph[["Family", "Level"]], left_on="fam1", right_on="Family", how='left')
predeats_ = predeats_.rename(columns={'Level': 'fam1_troph'})
predeats_ = pd.merge(predeats_, troph[["Family", "Level"]], left_on="fam2", right_on="Family", how='left')
predeats_ = predeats_.rename(columns={'Level': 'fam2_troph'})
predeats_ = predeats_.drop(columns=['Family_x', 'Family_y', 'interaction_type'])
predeats_ = predeats_.drop_duplicates()
exclude_rows = (predeats_["fam1"] == "Fragilariaceae") & (predeats_["fam2"] == "Aplysiidae") # the result of the two families is a database lookup mismatch
predeats_ = predeats_[~exclude_rows]

# identify the trophic level of the interaction. two interacting families can have different trophic levels. we want to assign the interaction to higher level.
grouped = predeats_.groupby(['fam1', 'fam2']).agg(
    count=('fam1', 'size'),                   # Count occurrences
    max_troph1=('fam1_troph', 'max'),          # Maximum of fam1_troph
    max_troph2=('fam2_troph', 'max')            # Maximum of fam2_troph
).reset_index()

# Compute the maximum between 'max_troph1' and 'max_troph2' for each group
grouped['troph'] = grouped[['max_troph1', 'max_troph2']].max(axis=1)
grouped_preds = grouped.drop(columns=['max_troph1', 'max_troph2'])

G_globi = nx.Graph()
for i in grouped_preds.index:
    edge_0 = grouped_preds.iloc[i]["fam1"]
    edge_1 = grouped_preds.iloc[i]["fam2"]
    edge_count = grouped_preds.iloc[i]["count"]
    trophic_level = grouped_preds.iloc[i]["troph"]
    G_globi.add_edge(edge_0, edge_1, weight=edge_count, troph=trophic_level)

# load the consensus network
cn = pd.read_csv(outdir+"KL77_cn.csv", delimiter=";", header=0, index_col=0)
G_cn = nx.from_pandas_adjacency(cn)

# load the weighted spiec easi network to get the same amount of edges
SE = load_spiec_easi(outdir+"KL-77_spieceasi_weighted")
se_thresh = 0.3365 # this empirical threshold results in 212 edges
se = SE.where(abs(SE) > se_thresh, 0.0)
se = se.where(se == 0, 1)
G_se = nx.from_pandas_adjacency(se)

# get the interactions which occur in the networks and save the interactions of the consensus network
cn_globi_ret = getInteractionDB(cn, predeats_)
se_globi_ret = getInteractionDB(se, predeats_)
cn_globi_ret.to_csv(outdir+"cn_globi_list.csv")

# get the trophic levels of the interactions and sum the occurences for each level
cn_troph_levels = []
for i in cn_globi_ret.index:
    cn_troph_levels.append(max(cn_globi_ret.iloc[i]["fam1_troph"], cn_globi_ret.iloc[i]["fam2_troph"]))

se_troph_levels = []
for i in se_globi_ret.index:
    se_troph_levels.append(max(se_globi_ret.iloc[i]["fam1_troph"], se_globi_ret.iloc[i]["fam2_troph"]))
   
cn_trophs = Counter(cn_troph_levels)
se_trophs = Counter(se_troph_levels)

## FAMILY COUNTS PER TROHIC LAYER
'''
troph["Family"].values
M = pd.read_csv(indir+"F_KL-77.csv")
new_columns = troph["Family"].values
test = troph.T
test.columns = new_columns
test = test.drop("Family")

matching_levels = test.reindex(columns=M.columns[1:]).fillna(0)
M.loc['Level'] = matching_levels.iloc[1]
M = pd.concat([M.loc[['Level']], M.loc[M.index != 'Level']])
M.to_csv("C:/Users/vdinkel/Desktop/famtroph.csv")
'''

def z_score(values, x):
    values = np.array(values)
    z = (x - np.mean(values)) / np.std(values)
    return round(z,2)

def dfs_longest_path(graph, current_node, visited_edges, path):
    longest_path = path[:]
    
    for neighbor in graph.neighbors(current_node):
        edge = frozenset([current_node, neighbor])  # Undirected edge
        if edge not in visited_edges:
            visited_edges.add(edge)  # Mark edge as visited
            candidate_path = dfs_longest_path(graph, neighbor, visited_edges, path + [neighbor])
            if len(candidate_path) > len(longest_path):
                longest_path = candidate_path
            visited_edges.remove(edge)  # Backtrack
    
    return longest_path

def empiricalP(null_values, T_obs):
    p_value = np.sum(null_values >= T_obs) / len(null_values)
    return round(p_value,2 )

# Find the longest path in the graph
def find_longest_chain(graph):
    longest_path = []
    for node in graph.nodes():
        visited_edges = set()  # Track visited edges
        candidate = dfs_longest_path(graph, node, visited_edges, [node])
        if len(candidate) > len(longest_path):
            longest_path = candidate
    return longest_path

def get_shortest_paths(G):
    lccs = [c for c in sorted(nx.connected_components(G), key=len, reverse=True)]
    shortest_paths = []
    longest_chains = []
    longest_chains_explicit = []
    diameters = []
    troph_distances = []
    lcc_sizes = []
    for lcc in lccs:
        G_sub = G.subgraph(lcc)
        shortest_paths.append(nx.average_shortest_path_length(G_sub))
        # Compute the longest chain
        longest_chain = find_longest_chain(G_sub)
        longest_chains_explicit.append(longest_chain)
        longest_chains.append(max([len(k) for k in longest_chain]))
        
        lcc_diameter = nx.diameter(G_sub)
        diameters.append(lcc_diameter)
        
        # trophic distance (maxtroph - mintroph)
        lcc_troph = troph[troph["Family"].isin(lcc)]["Level"].values
        troph_dist = max(lcc_troph) - min(lcc_troph)
        troph_distances.append(troph_dist)
        
        lcc_sizes.append(len(lcc))
        
    return troph_distances, diameters, shortest_paths, sorted(longest_chains, reverse=True), longest_chains_explicit, lcc_sizes

# get globi interactions for cn and se networks
G_cn, cn_globi_edge_overlap, cn_globi_troph_counts = mergeInteractions(G_cn, G_globi)
G_se, se_globi_edge_overlap, se_globi_troph_counts = mergeInteractions(G_se, G_globi)

# randomization for null model
rnd_edges = []
rnd_trophs = []

# get only the occuring interactions and the corresponding interaction graph
cn_interaction_edges = [edge for edge in G_cn.edges if G_cn.edges[edge[0], edge[1]]["weight"] > 0]
G_cn_interactions = G_cn.edge_subgraph(cn_interaction_edges)
se_interaction_edges = [edge for edge in G_se.edges if G_se.edges[edge[0], edge[1]]["weight"] > 0]
G_se_interactions = G_se.edge_subgraph(se_interaction_edges)

# get statistics of the interaction graphs
cn_troph_distances, cn_diamters, cn_shortest_paths, cn_longest_chains, cn_longest_chains_explicit, cn_lcc_sizes = get_shortest_paths(G_cn_interactions)
se_troph_distances, se_diamters, se_shortest_paths, se_longest_chains, se_longest_chains_explicit, se_lcc_sizes = get_shortest_paths(G_se_interactions)

# get interaction length = 2
def get_indirect_interactions(G, Globi_G):
    ind_edges = {}
    ind_interactions = {}
    for edge in G.edges:
        if edge not in Globi_G.edges():
            try:
                e0 = list(dict(Globi_G[edge[0]]).keys())
                e1 = list(dict(Globi_G[edge[1]]).keys())
                se01 = set(e0) & set(e1)
                if len(se01) > 0:
                    ind_edges[edge] = len(se01)
                    totsum = np.sum([Globi_G.edges[edge[0], k]["weight"] for k in list(se01)]) + np.sum([Globi_G.edges[edge[1], k]["weight"] for k in list(se01)])
                    ind_interactions[edge] = totsum
                    
                    e0_troph = troph[troph["Family"]== edge[0]]["Level"].values[0]
                    e1_troph = troph[troph["Family"]== edge[1]]["Level"].values[0]
                    #('Stephanodiscaceae', 'Blenniidae')
                    #if min(e0_troph, e1_troph)  == 1 and max(e0_troph , e1_troph) > 1:
                    #    import pdb; pdb.set_trace()
            except:
                pass
    return ind_edges, ind_interactions


cn_indirs, cn_ind_interactions = get_indirect_interactions(G_cn, G_globi) # menge der indirekten verbindungen zwischen einem paar, menge der interaktionen abgedeckt von den indirekten verbindungen
cn_indir_edges = len(cn_indirs.keys()) # wie viele kanten haben keine direkte verbindung in globi, dafür aber gemeinsame nachbarn? 81
cn_indir_edges_sum = np.sum([cn_indirs[k] for k in cn_indirs.keys()]) # was ist die summe dieser gemeinsamen nachbarn?
cn_indir_interactions_sum = np.sum([cn_ind_interactions[k] for k in cn_ind_interactions.keys()]) # wie viele interaktionen werden durch diese indirekten verbindungen abgedeckt?

se_indirs, se_ind_interactions = get_indirect_interactions(G_se, G_globi) # menge der indirekten verbindungen zwischen einem paar, menge der interaktionen abgedeckt von den indirekten verbindungen
se_indir_edges = len(se_indirs.keys()) # wie viele kanten haben keine direkte verbindung in globi, dafür aber gemeinsame nachbarn? 81
se_indir_edges_sum = np.sum([se_indirs[k] for k in se_indirs.keys()]) # was ist die summe dieser gemeinsamen nachbarn?
se_indir_interactions_sum = np.sum([se_ind_interactions[k] for k in se_ind_interactions.keys()]) # wie viele interaktionen werden durch diese indirekten verbindungen abgedeckt?

n_rnd_runs = 100 # randomization runs (10000)
rnd_interaction_edges = []
rnd_runs = []
rnd_indir_edges = []
rnd_indir_edges_sum = []
rnd_indir_interactions_sum = []

for i in range(0, n_rnd_runs):
    G_rnd = nx.from_pandas_adjacency(cn)
    nx.double_edge_swap(G_rnd, nswap=1000, max_tries=10000)
    
    G_rnd, edge_overlap, troph_counts  = mergeInteractions(G_rnd, G_globi) # new structure, since the Graph is now modified
    rnd_edges.append(edge_overlap)
    rnd_trophs.append(troph_counts)
    
    rnd_interaction_edges = [edge for edge in G_rnd.edges if G_rnd.edges[edge[0], edge[1]]["weight"] > 0]
    G_rnd_interactions = G_rnd.edge_subgraph(rnd_interaction_edges)
    #rnd_troph_distances, rnd_diameters, rnd_shortest_paths, rnd_longest_chains, rnd_longest_chains_explicit = get_shortest_paths(G_rnd_interactions)
    rnd_runs.append(get_shortest_paths(G_rnd_interactions))
    
    rnd_indir, rnd_ind_interactions = get_indirect_interactions(G_rnd, G_globi)
    rnd_indir_edges.append(len(rnd_indir.keys()))
    rnd_indir_edges_sum.append(np.sum([rnd_indir[k] for k in rnd_indir.keys()]))
    rnd_indir_interactions_sum.append(np.sum([rnd_ind_interactions[k] for k in rnd_ind_interactions.keys()]))
    
    if i % 500 == 0:
        print(f"Iteration {i}")

empiricalP(np.array(rnd_indir_edges), cn_indir_edges)# >0
empiricalP(np.array(rnd_indir_edges_sum), cn_indir_edges_sum) #= 0
empiricalP(np.array(rnd_indir_interactions_sum), cn_indir_interactions_sum) # = 0
# --> Verbundene Kanten im CN ergeben sich durch indirekte abhängigkeiten zu gemeinsamen nachbarn

# store randomization runs. comment the randomization out if it was done already and just load the stored files
pd.DataFrame(rnd_edges).to_csv(outdir+"/rnd_interaction_edges.csv")
pd.DataFrame(rnd_trophs).to_csv(outdir+"/rnd_interaction_trophs.csv")

# load randomization results
pd_rnd_edges = pd.read_csv(outdir+"/rnd_interaction_edges.csv")
pd_rnd_trophs = pd.read_csv(outdir+"/rnd_interaction_trophs.csv")

# cn_globi_ret, sn_globi_ret are the lists of Globi interactions
exportToGephi(G_cn, outdir+"cn_globi_gephi.gml", troph, cn_globi_ret, cn_indirs, G_man1)

# FIRST PLOT
fig, axs = plt.subplots(1, 3, sharey=False, tight_layout=True, figsize=(19, 6))

first_bins = 16
linestyle_cn = "-"
if cn_globi_edge_overlap == se_globi_edge_overlap:
    linestyle_cn = "--"
axs[0].hist(pd_rnd_edges["0"].values, color='lightgray', bins=first_bins)
axs[0].axvline(x=cn_globi_edge_overlap, color='red', linestyle=linestyle_cn, linewidth=2, label='Observed Value - Consensus') 
axs[0].axvline(x=se_globi_edge_overlap, color='green', linestyle='--', linewidth=2, label='Observed Value - SPIEC-EASI') 
axs[0].set_title("Total Overlap of Network Edges")
axs[0].set_xlabel("Overlapping Edges with Globi Network")
axs[0].set_ylabel("Frequency of Overlapping Edges/Taxa")

rnd_troph_sums = pd_rnd_trophs[["1", "2", "3", "4"]].T.sum().values
axs[1].hist(rnd_troph_sums, color='lightgray', bins=first_bins)
axs[1].axvline(x=sum(cn_globi_troph_counts.values()), color='red', linestyle='-', linewidth=2, label='Observed Value - Consensus')
axs[1].axvline(x=sum(se_globi_troph_counts.values()), color='green', linestyle='--', linewidth=2, label='Observed Value - SPIEC-EASI')
axs[1].set_title("Trophic Level 1-4 Taxa")
axs[1].set_xlabel("Overlapping Taxa with Globi Network")
axs[1].set_ylabel("Frequency of Overlapping Edges/Taxa")

i = 4
axs[2].hist(pd_rnd_trophs[str(i)].values, color='lightgray', bins=first_bins)
axs[2].axvline(x=cn_globi_troph_counts[str(i)], color='red', linestyle='-', linewidth=2, label='Observed Value - Consensus')
axs[2].axvline(x=se_globi_troph_counts[str(i)], color='green', linestyle='--', linewidth=2, label='Observed Value - SPIEC-EASI')
axs[2].set_title("Trophic Level "+str(i)+" Taxa")
axs[2].set_xlabel("Overlapping Taxa with Globi Network")
axs[2].set_ylabel("Frequency of Overlapping Edges/Taxa")

# FIRST - EDGE COUNTS
p_val_cn = empiricalP(pd_rnd_edges["0"].values, cn_globi_edge_overlap)
p_val_se = empiricalP(pd_rnd_edges["0"].values, se_globi_edge_overlap)
z_score_cn = z_score(pd_rnd_edges["0"].values, cn_globi_edge_overlap)
z_val_se = z_score(pd_rnd_edges["0"].values, se_globi_edge_overlap)
print("FIRST - EDGE COUNTS")
print("CN | z_score: ", z_score_cn, " p:", p_val_cn)
print("SE | z_score: ", z_val_se, " p:", p_val_se)

# FIRST - ALL TROPHIC (1-4) COUNTS
p_val_cn = empiricalP(rnd_troph_sums, sum(cn_globi_troph_counts.values()))
p_val_se = empiricalP(rnd_troph_sums, sum(se_globi_troph_counts.values()))
z_score_cn = z_score(rnd_troph_sums, sum(cn_globi_troph_counts.values()))
z_val_se = z_score(rnd_troph_sums, sum(se_globi_troph_counts.values()))
print("FIRST - ALL TROPHIC (1-4) COUNTS")
print("CN | z_score: ", z_score_cn, " p:", p_val_cn)
print("SE | z_score: ", z_val_se, " p:", p_val_se)

# FIRST - TROPHIC LEVEL 4 COUNTS
p_val_cn = empiricalP(pd_rnd_trophs[str(4)].values, cn_globi_troph_counts[str(4)])
p_val_se = empiricalP(pd_rnd_trophs[str(4)].values, se_globi_troph_counts[str(4)])
z_score_cn = z_score(pd_rnd_trophs[str(4)].values, cn_globi_troph_counts[str(4)])
z_val_se = z_score(pd_rnd_trophs[str(4)].values, se_globi_troph_counts[str(4)])
print("FIRST - ALL TROPHIC (1-4) COUNTS")
print("CN | z_score: ", z_score_cn, " p:", p_val_cn)
print("SE | z_score: ", z_val_se, " p:", p_val_se)

import matplotlib.lines as mlines
# Create custom legend entries
legend_handles = [
    mlines.Line2D([], [], color='lightgray', linewidth=6, label='Null Distribution'),
    mlines.Line2D([], [], color='red', linestyle='--', linewidth=1.5, label='Consensus Network'),
    mlines.Line2D([], [], color='green', linestyle='--', linewidth=1.5, label='SPIEC-EASI')
]

# Add a single legend outside the subplots
fig.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=3)
plt.savefig(plotdir+"S4_troph_overlap_A.png")
plt.show()

##### SECOND PLOT - this was removed
'''
second_bins = 16
fig, axs = plt.subplots(1, 3, sharey=False, tight_layout=True, figsize=(19, 6))

axs[0].hist([np.mean(rnd_runs[k][5]) for k in range(0,n_rnd_runs)], color='lightgray', bins=second_bins)
axs[0].axvline(x=np.mean(cn_lcc_sizes), color='red', linestyle="-", linewidth=2, label='Observed Value - Consensus') 
axs[0].axvline(x=np.mean(se_lcc_sizes), color='green', linestyle='--', linewidth=2, label='Observed Value - SPIEC-EASI') 
axs[0].set_title("Mean LCC Sizes (Interactions)")
axs[0].set_xlabel("Mean LCC Size")
axs[0].set_ylabel("Frequency of Sizes")

z_val_1 = z_score([np.mean(rnd_runs[k][5]) for k in range(0,n_rnd_runs)], np.mean(cn_lcc_sizes))
p_val_1 = empiricalP([np.mean(rnd_runs[k][5]) for k in range(0,n_rnd_runs)], np.mean(cn_lcc_sizes))

axs[1].hist([np.mean(rnd_runs[k][1]) for k in range(0,n_rnd_runs)], color='lightgray', bins=second_bins)
axs[1].axvline(x=np.mean(cn_diamters), color='red', linestyle='-', linewidth=2, label='Observed Value - Consensus')
axs[1].axvline(x=np.mean(se_diamters), color='green', linestyle='--', linewidth=2, label='Observed Value - SPIEC-EASI')
axs[1].set_title("Diameters")
axs[1].set_xlabel("Mean Diameter")
axs[1].set_ylabel("Frequency of Diameters")

z_val_2 = z_score([np.mean(rnd_runs[k][1]) for k in range(0,n_rnd_runs)], np.mean(cn_diamters))
axs[1].text(
    0.97, 0.99, f'z-score CN = {z_val_2}', color='red', ha='right', va='top', transform=axs[1].transAxes, fontsize=10
)
p_val_2 = empiricalP([np.mean(rnd_runs[k][1]) for k in range(0,n_rnd_runs)], np.mean(cn_diamters))
axs[1].text(
    0.97, 0.95, f'p-value CN = {p_val_2}', color='red', ha='right', va='top', transform=axs[1].transAxes, fontsize=10
)

axs[2].hist([np.mean(rnd_runs[k][2]) for k in range(0,n_rnd_runs)], color='lightgray', bins=second_bins)
axs[2].axvline(x=np.mean(cn_shortest_paths), color='red', linestyle='-', linewidth=2, label='Observed Value - Consensus')
axs[2].axvline(x=np.mean(se_shortest_paths), color='green', linestyle='--', linewidth=2, label='Observed Value - SPIEC-EASI')
axs[2].set_title("Shortest Paths")
axs[2].set_xlabel("Mean Shortest Paths")
axs[2].set_ylabel("Frequency of Shortest Paths")

z_val_3 = z_score([np.mean(rnd_runs[k][2]) for k in range(0,n_rnd_runs)], np.mean(cn_shortest_paths))
axs[2].text(
    0.97, 0.99, f'z-score CN = {z_val_3}', color='red', ha='right', va='top', transform=axs[2].transAxes, fontsize=10
)
p_val_3 = empiricalP([np.mean(rnd_runs[k][2]) for k in range(0,n_rnd_runs)], np.mean(cn_shortest_paths))
axs[2].text(
    0.97, 0.95, f'p-value CN = {p_val_3}', color='red', ha='right', va='top', transform=axs[2].transAxes, fontsize=10
)

for i in range(0,3):
    axs[i].spines['top'].set_visible(False)    # Remove the top border
    axs[i].spines['right'].set_visible(False)  # Remove the right border
    axs[i].spines['left'].set_visible(False)   # Remove the left border (optional)
    axs[i].spines['bottom'].set_visible(False) # Remove the bottom border (optional)

import matplotlib.lines as mlines
# Create custom legend entries
legend_handles = [
    mlines.Line2D([], [], color='lightgray', linewidth=6, label='Null Distribution'),
    mlines.Line2D([], [], color='red', linestyle='--', linewidth=1.5, label='Consensus Network'),
    mlines.Line2D([], [], color='green', linestyle='--', linewidth=1.5, label='SPIEC-EASI')
]

# Add a single legend outside the subplots
fig.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=3)
plt.show()
'''


##### THIRD PLOT
third_bins = 16
fig, axs = plt.subplots(1, 3, sharey=False, tight_layout=True, figsize=(19, 6))

axs[0].hist(rnd_indir_edges, color='lightgray', bins=third_bins)
axs[0].axvline(x=cn_indir_edges, color='red', linestyle="-", linewidth=2, label='Observed Value - Consensus') 
axs[0].axvline(x=se_indir_edges, color='green', linestyle='--', linewidth=2, label='Observed Value - SPIEC-EASI') 
axs[0].set_title("Indirect Edges")
axs[0].set_xlabel("Indirect Edges")
axs[0].set_ylabel("Frequency")

axs[1].hist(rnd_indir_edges_sum, color='lightgray', bins=third_bins)
axs[1].axvline(x=cn_indir_edges_sum, color='red', linestyle='-', linewidth=2, label='Observed Value - Consensus')
axs[1].axvline(x=se_indir_edges_sum, color='green', linestyle='--', linewidth=2, label='Observed Value - SPIEC-EASI')
axs[1].set_title("Number of Indirect Families")
axs[1].set_xlabel("Number of Indirect Families")
axs[1].set_ylabel("Frequency")

axs[2].hist(rnd_indir_interactions_sum, color='lightgray', bins=third_bins)
axs[2].axvline(x=cn_indir_interactions_sum, color='red', linestyle='-', linewidth=2, label='Observed Value - Consensus')
axs[2].axvline(x=se_indir_interactions_sum, color='green', linestyle='--', linewidth=2, label='Observed Value - SPIEC-EASI')
axs[2].set_title("Sum of Indirect Interactions")
axs[2].set_xlabel("Sum of Indirect Interactions")
axs[2].set_ylabel("Frequency")

import matplotlib.lines as mlines
# Create custom legend entries
legend_handles = [
    mlines.Line2D([], [], color='lightgray', linewidth=6, label='Null Distribution'),
    mlines.Line2D([], [], color='red', linestyle='--', linewidth=1.5, label='Consensus Network'),
    mlines.Line2D([], [], color='green', linestyle='--', linewidth=1.5, label='SPIEC-EASI')
]

# Add a single legend outside the subplots
fig.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=3)
plt.savefig(plotdir+"S4_troph_overlap_B.png")
plt.show()

# THIRD - INDIRECT EDGE COUNTS
p_val_cn = empiricalP(np.array(rnd_indir_edges), cn_indir_edges)
p_val_se = empiricalP(np.array(rnd_indir_edges), se_indir_edges)
z_score_cn = z_score(rnd_indir_edges, cn_indir_edges)
z_val_se = z_score(rnd_indir_edges, se_indir_edges)
print("THIRD - INDIRECT EDGE COUNTS")
print("CN | z_score: ", z_score_cn, " p:", p_val_cn)
print("SE | z_score: ", z_val_se, " p:", p_val_se)

# THIRD - INDIRECT FAMILY COUNTS
p_val_cn = empiricalP(np.array(rnd_indir_edges_sum), cn_indir_edges_sum)
p_val_se = empiricalP(np.array(rnd_indir_edges_sum), se_indir_edges_sum)
z_score_cn = z_score(rnd_indir_edges_sum, cn_indir_edges_sum)
z_val_se = z_score(rnd_indir_edges_sum, se_indir_edges_sum)
print("THIRD - INDIRECT EDGE COUNTS")
print("CN | z_score: ", z_score_cn, " p:", p_val_cn)
print("SE | z_score: ", z_val_se, " p:", p_val_se)

# THIRD - SUM OF INDIRECT TAXA
p_val_cn = empiricalP(np.array(rnd_indir_interactions_sum), cn_indir_interactions_sum)
p_val_se = empiricalP(np.array(rnd_indir_interactions_sum), se_indir_interactions_sum)
z_score_cn = z_score(rnd_indir_interactions_sum, cn_indir_interactions_sum)
z_val_se = z_score(rnd_indir_interactions_sum, se_indir_interactions_sum)
print("THIRD - INDIRECT EDGE COUNTS")
print("CN | z_score: ", z_score_cn, " p:", p_val_cn)
print("SE | z_score: ", z_val_se, " p:", p_val_se)