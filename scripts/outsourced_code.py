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

# SANKEY DIAGRAM PLOT
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