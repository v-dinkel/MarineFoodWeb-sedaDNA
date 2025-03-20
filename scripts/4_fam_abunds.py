# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 01:54:14 2023

@author: vdinkel
"""


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import json


def read_config(filename):
    #read the config file of the project
    config = {}
    with open(filename, 'r') as file:
        for line in file:
            if line.strip() and not line.startswith("#"):  # Ignore empty lines or comments
                key, value = line.strip().split(":", 1)
                config[key.strip()] = value.strip()
    return config

def loadFullKrakenReport(krakenFileDir, hasHeader):
    print("--loadFullKrakenReport")
    if hasHeader:
        headerID = 0
    else:
        headerID = None
    df = pd.read_table(krakenFileDir, sep=",", header=headerID)
    return df


workdir = read_config("../config.txt")["workdir"]
outdir = workdir+"output/"
suppdir = workdir+"supplementary_information/"
abunds = workdir+"input/F_KL-77.csv"
rel_abunds =  workdir+"output/F_KL-77_rel.csv"
fam_troph = suppdir+"families_trophic_functions.csv"
libs =  workdir+"input/KL-77_nt0.2.csv"

# get positive / negative ngrip correlations
scorrs = pd.read_table(outdir+"KL-77_temp_spearman_corrs.csv", sep=";", header=1)
scorrs_p = pd.read_table(outdir+"KL-77_temp_spearman_corrs_p.csv", sep=";", header=1)
sc_df = scorrs.T.iloc[1:]
tpos2 = list(sc_df[sc_df > 0.3].dropna().index)
tneg2 = list(sc_df[sc_df < -0.3].dropna().index)

with open(outdir+"troph_level_dict.json", 'r') as json_file:
    trDict = json.load(json_file)

#df_krakenReport = loadFullKrakenReport(libs, hasHeader = True)
df_abunds = pd.read_table(abunds, sep=",", header=0)
df_rel_abunds = pd.read_table(rel_abunds, sep=";")
df_fam_troph =  pd.read_table(fam_troph, sep=",")

# statistics of classified read counts etc. requires full kraken report
'''
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
'''

# STRATIGRAM
ages = df_abunds["Unnamed: 0"].values
ngrip = pd.read_table(workdir+"supplementary_information/families_ngrip_corrs.csv", sep=",", header=0)
keystones = pd.read_table(workdir+"supplementary_information/cn_centralities.csv", sep=";", header=0)

#trColors = ["#6c655a", "#363a15", "#485c6fa", "#8a5927"]
trColors = ["#84C03A", "#415b14", "#273f60", "#855f06"]

allgrps = []
nplots = 0
for trlvl in range(1,5):
    
    occ_labels = list(keystones["Unnamed: 0"].values)
    intersect = list(set(trDict[str(trlvl)]) & set(occ_labels))

    merged_labels = list(df_rel_abunds[intersect].sum().sort_values(ascending=False).index[:15])
    sm_grp_abunds = df_rel_abunds[merged_labels]
    allgrps.append(sm_grp_abunds)
    nplots += len(sm_grp_abunds.columns)

fig, axs = plt.subplots(nplots,1, figsize=(12,22), sharex=True, sharey=False)

i = 0
k = 0

for gr in allgrps:
    for col in gr.columns:

        scolor = "#000000"
        if col in ngrip["Unnamed: 0"].values:
            if ngrip[ngrip["Unnamed: 0"] == col]["ngrip_r"].values < 0:
                scolor = "#0077BB"
            else:
                scolor = "#CC3311"
        
        axs[i].spines['top'].set_visible(False)
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
        
        yticks = axs[i].yaxis.get_major_ticks()
        for ytick in yticks:
            ytick.label.set_label("")
            ytick.label1.set_label("")
            ytick.label2.set_label("")
        
        yy = True
        for n, label in enumerate(axs[i].yaxis.get_ticklabels()):
            if n == 0:
                label.set_visible(False)
        
        i += 1
    k += 1

plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.1, wspace=0.0, hspace=0)

plt.savefig(workdir+'plots/family_abundances.png', format='png', dpi=1200)
plt.show()