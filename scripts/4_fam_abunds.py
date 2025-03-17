# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 01:54:14 2023

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

def loadFullKrakenReport(krakenFileDir, hasHeader):
    print("--loadFullKrakenReport")
    if hasHeader:
        headerID = 0
    else:
        headerID = None
    df = pd.read_table(krakenFileDir, sep=",", header=headerID)
    return df

tpos = ['Acipenseridae',
 'Acroporidae',
 'Balaenopteridae',
 'Batrachoididae',
 'Bolbocoleonaceae',
 'Chaetocerotaceae',
 'Cionidae',
 'Clupeidae',
 'Collosphaeridae',
 'Crustomastigaceae',
 'Delphinidae',
 'Echeneidae',
 'Fundulidae',
 'Gasterosteidae',
 'Geminigeraceae',
 'Gonyaulacaceae',
 'Harpacticidae',
 'Kallymeniaceae',
 'Lateolabracidae',
 'Mustelidae',
 'Mytilidae',
 'Pelagiidae',
 'Penaeidae',
 'Phocoenidae',
 'Salmonidae',
 'Solieriaceae',
 'Sparidae',
 'Triparmaceae',
 'Ulmaridae']
tneg = ['Acrochaetiaceae',
 'Aphanizomenonaceae',
 'Areschougiaceae',
 'Attheyaceae',
 'Bacillariaceae',
 'Calotrichaceae',
 'Chlamydomonadaceae',
 'Cirratulidae',
 'Cottidae',
 'Cyprinidae',
 'Euplotidae',
 'Gadidae',
 'Histionidae',
 'Monodopsidaceae',
 'Parastacidae',
 'Phaeocystaceae',
 'Picocystaceae',
 'Polynoidae',
 'Priapulidae',
 'Pyuridae',
 'Suessiaceae',
 'Zygnemataceae']
trDict = {'1': ['Acanthocerataceae',
  'Acaryochloridaceae',
  'Acinetosporaceae',
  'Acrochaetiaceae',
  'Agaraceae',
  'Alariaceae',
  'Amoebophryaceae',
  'Amphipleuraceae',
  'Anaulaceae',
  'Anomoeoneidaceae',
  'Aphanizomenonaceae',
  'Aphanothecaceae',
  'Areschougiaceae',
  'Attheyaceae',
  'Bacillariaceae',
  'Bangiaceae',
  'Bathycoccaceae',
  'Batrachospermaceae',
  'Berkeleyaceae',
  'Biddulphiaceae',
  'Bolbocoleonaceae',
  'Bonnemaisoniaceae',
  'Bracteacoccaceae',
  'Bryopsidaceae',
  'Calotrichaceae',
  'Ceramiaceae',
  'Chaetocerotaceae',
  'Chamaesiphonaceae',
  'Characeae',
  'Chattonellaceae',
  'Chlamydomonadaceae',
  'Chlorellaceae',
  'Chlorobiaceae',
  'Chlorocystidaceae',
  'Chloroflexaceae',
  'Chloropicaceae',
  'Chromulinaceae',
  'Chroococcaceae',
  'Chroococcidiopsidaceae',
  'Chrysochromulinaceae',
  'Chrysolepidomonadaceae',
  'Closteriaceae',
  'Codiaceae',
  'Coleofasciculaceae',
  'Corallinaceae',
  'Coscinodiscaceae',
  'Crustomastigaceae',
  'Cyanidiaceae',
  'Cyanophoraceae',
  'Cyanothecaceae',
  'Cymatosiraceae',
  'Dasyaceae',
  'Delesseriaceae',
  'Derbesiaceae',
  'Dermocarpellaceae',
  'Desmidiaceae',
  'Dictyotaceae',
  'Dunaliellaceae',
  'Ectocarpaceae',
  'Entomoneidaceae',
  'Eunotiaceae',
  'Eustigmataceae',
  'Fragilariaceae',
  'Fucaceae',
  'Gelidiaceae',
  'Gigartinaceae',
  'Gloeobacteraceae',
  'Gloeomargaritaceae',
  'Gomontiellaceae',
  'Gomphonemataceae',
  'Gonyaulacaceae',
  'Gracilariaceae',
  'Gymnodiniaceae',
  'Haematococcaceae',
  'Halimedaceae',
  'Halymeniaceae',
  'Hapalidiaceae',
  'Hapalosiphonaceae',
  'Heliobacteriaceae',
  'Hemiaulaceae',
  'Hemidiscaceae',
  'Heterocapsaceae',
  'Hildenbrandiaceae',
  'Hydrodictyaceae',
  'Hyellaceae',
  'Hymenomonadaceae',
  'Hypneaceae',
  'Kallymeniaceae',
  'Kareniaceae',
  'Klebsormidiaceae',
  'Koliellaceae',
  'Kryptoperidiniaceae',
  'Laminariaceae',
  'Leptocylindraceae',
  'Leptolyngbyaceae',
  'Liagoraceae',
  'Licmophoraceae',
  'Lithodesmiaceae',
  'Mallomonadaceae',
  'Mamiellaceae',
  'Merismopediaceae',
  'Mesotaeniaceae',
  'Microcoleaceae',
  'Monodopsidaceae',
  'Naviculaceae',
  'Nephroselmidaceae',
  'Noelaerhabdaceae',
  'Nostocaceae',
  'Oedogoniaceae',
  'Oltmannsiellopsidaceae',
  'Oocystaceae',
  'Oscillatoriaceae',
  'Ostreobiaceae',
  'Palmariaceae',
  'Pavlovaceae',
  'Pfiesteriaceae',
  'Phaeocystaceae',
  'Phaeodactylaceae',
  'Phyllophoraceae',
  'Picocystaceae',
  'Pinguiochrysidaceae',
  'Plagiogrammaceae',
  'Pleurastraceae',
  'Polyphysaceae',
  'Porphyridiaceae',
  'Prasinococcaceae',
  'Prasiolaceae',
  'Prochloraceae',
  'Prochlorotrichaceae',
  'Prorocentraceae',
  'Protoperidiniaceae',
  'Pseudanabaenaceae',
  'Pycnococcaceae',
  'Pyramimonadaceae',
  'Rhizosoleniaceae',
  'Rhodogorgonaceae',
  'Rhodomelaceae',
  'Rhodymeniaceae',
  'Rivulariaceae',
  'Roseiflexaceae',
  'Sargassaceae',
  'Scenedesmaceae',
  'Scytonemataceae',
  'Scytosiphonaceae',
  'Selenastraceae',
  'Skeletonemataceae',
  'Solieriaceae',
  'Spyridiaceae',
  'Staurosiraceae',
  'Stephanodiscaceae',
  'Suessiaceae',
  'Syndiniaceae',
  'Synechococcaceae',
  'Thalassiosiraceae',
  'Thraustochytriaceae',
  'Toxariaceae',
  'Trebouxiaceae',
  'Triceratiaceae',
  'Triparmaceae',
  'Ulnariaceae',
  'Ulvaceae',
  'Wrangeliaceae',
  'Zosteraceae',
  'Zygnemataceae'],
 '2': ['Acroporidae',
  'Actiniidae',
  'Acytosteliaceae',
  'Aiptasiidae',
  'Apusomonadidae',
  'Arcidae',
  'Calanidae',
  'Campanulariidae',
  'Cavenderiaceae',
  'Chroomonadaceae',
  'Cionidae',
  'Cirratulidae',
  'Collosphaeridae',
  'Cyaneidae',
  'Diplonemidae',
  'Edwardsiidae',
  'Eirenidae',
  'Eucalanidae',
  'Euglenaceae',
  'Euphausiidae',
  'Euplotidae',
  'Fonticulaceae',
  'Geminigeraceae',
  'Globorotaliidae',
  'Gromiidae',
  'Hemiselmidaceae',
  'Hexamitidae',
  'Histionidae',
  'Hydridae',
  'Hypotrichomonadidae',
  'Jakobidae',
  'Lingulidae',
  'Merulinidae',
  'Metopidae',
  'Metridinidae',
  'Mytilidae',
  'Nephtheidae',
  'Niphatidae',
  'Oikopleuridae',
  'Ostreidae',
  'Oxytrichidae',
  'Parameciidae',
  'Paulinellidae',
  'Pectinidae',
  'Pelagiidae',
  'Perkinsidae',
  'Phacaceae',
  'Pocilloporidae',
  'Priapulidae',
  'Pteriidae',
  'Pyuridae',
  'Rotaliidae',
  'Saccosporidae',
  'Salpidae',
  'Salpingoecidae',
  'Sphaerozoidae',
  'Sticholonchidae',
  'Styelidae',
  'Temoridae',
  'Terebellidae',
  'Tetrahymenidae',
  'Thaumatomastigidae',
  'Ulmaridae',
  'Vahlkampfiidae'],
 '3': ['Acanthasteridae',
  'Acipenseridae',
  'Aegisthidae',
  'Alepocephalidae',
  'Ammodytidae',
  'Anarhichadidae',
  'Anoplopomatidae',
  'Aplysiidae',
  'Apogonidae',
  'Asteriidae',
  'Bathylagidae',
  'Batrachoididae',
  'Blenniidae',
  'Ceratiidae',
  'Clupeidae',
  'Cottidae',
  'Cyclopteridae',
  'Cynoglossidae',
  'Cyprinidae',
  'Daphniidae',
  'Dugesiidae',
  'Echeneidae',
  'Gadidae',
  'Gasterosteidae',
  'Gecarcinidae',
  'Gobiidae',
  'Gonostomatidae',
  'Haliotidae',
  'Harpacticidae',
  'Hyalellidae',
  'Lateolabracidae',
  'Limulidae',
  'Liparidae',
  'Loliginidae',
  'Lottiidae',
  'Lysianassidae',
  'Microstomatidae',
  'Myctophidae',
  'Myxinidae',
  'Octopodidae',
  'Oreosomatidae',
  'Osmeridae',
  'Oxystominidae',
  'Paralichthyidae',
  'Parastacidae',
  'Penaeidae',
  'Petromyzontidae',
  'Pleuronectidae',
  'Polynoidae',
  'Portunidae',
  'Provannidae',
  'Rajidae',
  'Rhincodontidae',
  'Salmonidae',
  'Sebastidae',
  'Sparidae',
  'Strongylocentrotidae',
  'Syngnathidae',
  'Tetragonicipitidae',
  'Zoarcidae'],
 '4': ['Balaenopteridae',
  'Delphinidae',
  'Monodontidae',
  'Mustelidae',
  'Odobenidae',
  'Otariidae',
  'Phocidae',
  'Phocoenidae',
  'Physeteridae']}

workdir = 'C:/Users/vdinkel/Desktop/Manuscript/'
abunds = workdir+"input/F_KL-77.csv"
rel_abunds =  workdir+"output/F_KL-77_rel.csv"
fam_troph = workdir+"supplementary_information/families_trophic_functions.csv"
libs =  workdir+"input/KL-77_nt0.2.csv"

df_krakenReport = loadFullKrakenReport(libs, hasHeader = True)
df_abunds = loadFullKrakenReport(abunds, hasHeader = True)
df_rel_abunds = pd.read_table(rel_abunds, sep=";")
df_fam_troph =  pd.read_table(fam_troph, sep=",")

C_sum_df = df_krakenReport[df_krakenReport["Rank"] == "R"][["CladeCount", "lib_id"]].groupby(by=["lib_id"]).sum()
U_sum_df = df_krakenReport[df_krakenReport["Rank"] == "U"][["CladeCount", "lib_id"]].groupby(by=["lib_id"]).sum()

famNames = list(df_abunds.columns[1:])
# sum of all families
F_all_sum = df_krakenReport[df_krakenReport["Rank"] == "F"][["CladeCount", "lib_id"]].groupby(by=["lib_id"]).sum()

# sum of only occuring families
F_occ_sum = df_krakenReport[df_krakenReport["Name"].isin(famNames)][["Name","CladeCount","lib_id"]].groupby("lib_id").sum()

allS = 0
trFams = []
trOccs = []
for i in range(1,5):
    tr_level_occs = df_fam_troph[df_fam_troph["trophic_level"] == i]["family"].values
    tr_sum = df_rel_abunds[tr_level_occs].sum().sum()
    allS += tr_sum
    trOccs.append(tr_sum)
    trFams.append(len(tr_level_occs))
# np.array(trOccs)/42


# compare => abundances are the same!
sorted(df_abunds[df_abunds.columns[1:]].sum(axis=1))
sorted(F_occ_sum["CladeCount"].values)

# make new (final) DF:
# lib_id, lib_size, unclassified, classified, families, filtered families
# read numbers = report.html value "passed filter" divided by 2
newDF = df_krakenReport[['lib_id', "age"]].groupby(by=["lib_id"]).mean()
newDF.columns = ["Age"]
newDF = newDF.assign(TotalReads=(C_sum_df + U_sum_df)["CladeCount"])
newDF = newDF.assign(Unclassified=U_sum_df["CladeCount"])
newDF = newDF.assign(Classified=C_sum_df["CladeCount"])
newDF = newDF.assign(Families=F_all_sum["CladeCount"])
newDF = newDF.assign(FilteredFamilies=F_occ_sum["CladeCount"])
newDF.sort_values("Age")
#newDF.to_csv('C:/Users/vdinkel/Desktop/Data/LOGIC_realnets/output/statistics/all_read_statistics.csv', index=True)  

# STRATIGRAM
#gephi_export_file = "C:/Users/vdinkel/Desktop/Data/LOGIC_realnets/output/statistics/gephi_data_export.csv"
#gephi_data = pd.read_table(gephi_export_file, sep=";")
ages = df_abunds["Unnamed: 0"].values
ngrip = pd.read_table(workdir+"supplementary_information/families_ngrip_corrs.csv", sep=",", header=0)
keystones = pd.read_table(workdir+"supplementary_information/cn_centralities.csv", sep=";", header=0)

#trColors = ["#6c655a", "#363a15", "#485c6fa", "#8a5927"]
trColors = ["#84C03A", "#415b14", "#273f60", "#855f06"]

allgrps = []
nplots = 0
for trlvl in range(1,5):
    #grp = gephi_data[gephi_data['trophic_level']==trlvl]["Label"].values
    #sgrp = gephi_data[(gephi_data['trophic_level']==trlvl) & (gephi_data['scorr']!=0)]["Label"].values
    #sm_grp_gephi = gephi_data[gephi_data["Label"].isin(grp)].sort_values(by="pageranks", ascending=False)[:10]
    
    # merge top 10 pageranks with pos/neg correlations, which were not included
    #merged_labels = list(set(list(sm_grp_gephi["Label"].values) + list(sgrp)))
    
    occ_labels = list(keystones["Unnamed: 0"].values)
    intersect = list(set(trDict[str(trlvl)]) & set(occ_labels))

    merged_labels = list(df_rel_abunds[intersect].sum().sort_values(ascending=False).index[:15])
    sm_grp_abunds = df_rel_abunds[merged_labels]
    allgrps.append(sm_grp_abunds)
    nplots += len(sm_grp_abunds.columns)

#nplots += 1
fig, axs = plt.subplots(nplots,1, figsize=(12,22), sharex=True, sharey=False)


#axs[0].plot(ngrip["age"].values, ngrip["temp"].values)
#axs[0].set_title("NGRIP", x=-0.04, y=-0.1, horizontalalignment="right")

i = 0
k = 0
#gephi_data = pd.read_table(gephi_export_file, sep=";")
for gr in allgrps:
    for col in gr.columns:
        #scorr = gephi_data[gephi_data["Label"]==col]["scorr"].values[0]
        
        scolor = "#000000"
        if col in ngrip["Unnamed: 0"].values:
            if ngrip[ngrip["Unnamed: 0"] == col]["ngrip_r"].values < 0:
                scolor = "#0077BB"
            else:
                scolor = "#CC3311"
        
        axs[i].spines['top'].set_visible(False)
        #axs[i].spines['left'].set_visible(False)
        axs[i].spines['bottom'].set_visible(False)
            
        axs[i].plot(ages, gr[col].values, color=trColors[k])
        axs[i].fill_between(ages, gr[col].values, color=trColors[k], alpha=.5)
        ff = 0
        
        for age in ages:
            axs[i].plot([age, age], [0, gr[col].values[ff]], color=trColors[k])
            ff+=1
        axs[i].set_title(col, color=scolor, x=-0.04, y=-0.1, horizontalalignment="right")
        axs[i].set_xticklabels(axs[i].get_xticks(), fontsize=12)
        axs[i].xaxis.set_major_formatter(plt.FormatStrFormatter("%.0f"))
        axs[i].yaxis.tick_right()
        
        labels = []
        for lab in axs[i].get_yticks():
            if lab == 0:
                labels.append("")
            else:
                labels.append(str(lab))

        axs[i].set_yticklabels(labels, fontsize=8)
        #axs[i].yaxis.set_major_formatter(plt.FormatStrFormatter("%.3f"))
        
        yticks = axs[i].yaxis.get_major_ticks()
        for ytick in yticks:
            ytick.label.set_label("")
            ytick.label1.set_label("")
            ytick.label2.set_label("")
        #import pdb; pdb.set_trace()
        
        #axs[i].yaxis.set_major_locator(plt.MaxNLocator(1))
        
        yy = True
        for n, label in enumerate(axs[i].yaxis.get_ticklabels()):
            #if n % every_nth != 0:
            #import pdb; pdb.set_trace()
            if n == 0:
                label.set_visible(False)
        
        i += 1
    k += 1
#plt.locator_params(axis='y', nbins=1)
plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.1, wspace=0.0, hspace=0)

plt.savefig(workdir+'plots/family_abundances.png', format='png', dpi=1200)
plt.show()
import pdb; pdb.set_trace()
# plt.plot(df_rel_abunds[trDict["1"]].sum(axis=1))

def getModuleMTL(module):
    dist = []
    dist_abs = []
    for i in range(1,5):
        module_fams = gephi_data[(gephi_data["module"]==module) & (gephi_data["trophic_level"]==i)]["Label"].values
        abund = df_rel_abunds[module_fams].sum().sum()
        dist.append(abund)
        dist_abs.append(len(module_fams))
    return dist, dist_abs

for k in range(0,4):
    mtl, mtl_abs = getModuleMTL(k)
    print("NORM ", k," = ", np.sum((mtl / np.sum(mtl)) * np.array([1,2,3,4])))
    print("ABS ", k, " = ", np.sum(mtl_abs / np.sum(mtl_abs) * np.array([1,2,3,4])))
    import pdb; pdb.set_trace()

#---------------------------------------------------
import plotly.graph_objects as go
source = []
target = []
value = []
link_cols = []
trophic_levels_f = "C:/Users/vdinkel/Desktop/Data/LOGIC_realnets/input/trophic_levels.csv"
trophic_levels = pd.read_table(trophic_levels_f, sep=",", header = None)

# WORKFLOW
labels = ["Taxonomic abundance", "Inference Methods Header", "List of Inference Methods", "Spiec-Easi", "Spearman", "ecoCopula","Esabo", "Sparcc", "CCREPE", "Propr", "Box_TN", "Box_TTN", "Box_Filter", "Box_CN", "Box_Out"]
source.append(0)
target.append(1)
value.append(8)

source.append(1)
target.append(2)
value.append(8)

source.append(2)
target.append(3)
value.append(2)

source.append(2)
target.append(4)
value.append(1)

source.append(2)
target.append(5)
value.append(1)

source.append(2)
target.append(6)
value.append(1)

source.append(2)
target.append(7)
value.append(1)

source.append(2)
target.append(8)
value.append(1)

source.append(2)
target.append(9)
value.append(1)

# Methods to first Box
source.append(3)
target.append(10)
value.append(1)

source.append(4)
target.append(10)
value.append(1)

source.append(5)
target.append(10)
value.append(1)

source.append(6)
target.append(10)
value.append(1)

source.append(7)
target.append(10)
value.append(1)

source.append(8)
target.append(10)
value.append(1)

source.append(9)
target.append(10)
value.append(1)

# Box to Box
source.append(10)
target.append(11)
value.append(7)

# Box to Box
source.append(11)
target.append(12)
value.append(7)

# Box to Box
source.append(12)
target.append(13)
value.append(8)

# SE to Filter
source.append(3)
target.append(12)
value.append(1)

# Box to Box
source.append(13)
target.append(14)
value.append(8)


'''
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
'''
'''
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
'''
node_colors = []
# SANKEY DIAGRAM EXPORT
fig = go.Figure(data=[go.Sankey(
    orientation = "h",
    node = dict(
      pad = 20,
      thickness = 40,
      #y= [0.5 for k in range(0,len(source))],      #x= [0.2 for i in range(0,4)]+[0.2 for i in range(0,4)]+[0.5 for k in range(8,len(source))],
      line = dict(color = "black", width = 1.1),
      #label = labels,
      color = node_colors
    ),
    link = dict(
      source = source, # indices correspond to labels, eg A1, A2, A1, B1, ...
      target = target,
      value = value,
      #color = link_cols
  ))])

fig.update_layout(
    hovermode = 'x',
    title="<b>Test",
    font=dict(size = 12, color = 'black'),
    plot_bgcolor='white',
    paper_bgcolor='white',
    )

fig.write_html("C:/Users/vdinkel/Desktop/Data/LOGIC_realnets/output/statistics/sankey-diagram-workflow.html")