# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 10:15:47 2025

@author: vdinkel
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

def read_config(filename):
    #read the config file of the project
    config = {}
    with open(filename, 'r') as file:
        for line in file:
            if line.strip() and not line.startswith("#"):  # Ignore empty lines or comments
                key, value = line.strip().split(":", 1)
                config[key.strip()] = value.strip()
    return config

workdir = read_config("../config.txt")["workdir"]

# NGRIP from manuscript
ngrip = pd.read_table(workdir+"output/KL-77_ngrip.csv", sep=";", header=0)

# Full NGRIP from source
ngrip_full = pd.read_table(workdir+"input/ngrip.csv")
ngrip_full.columns = ["Age", "Temp"]

# KL77 IP25 - plot centered IP25 values
kl77_ip25 = pd.read_csv(workdir+"input/IP25_77KL.csv", sep=";")
kl77_ip25_ngrip = pd.merge_asof(kl77_ip25, ngrip_full, on='Age', direction='nearest') # get nearest ngrip values to the timestamps of the ip25 values
kl77_ip25_ngrip.columns = ["Age", "IP25", "NGRIP"]

kl77_ip25 = kl77_ip25_ngrip["IP25"].values
kl77_ip25_centered = (kl77_ip25 - np.mean(kl77_ip25)) / np.std(kl77_ip25, ddof=0)
kl77_ngrip = kl77_ip25_ngrip["NGRIP"].values
kl77_ngrip_centered = (kl77_ngrip - np.mean(kl77_ngrip)) / np.std(kl77_ngrip, ddof=0)

plt.figure(figsize=(8, 6), dpi=300)
plt.plot(kl77_ip25_ngrip["Age"]*1000, kl77_ip25_centered, label="KL77 IP$_{25}$", color="#1f77b4")
plt.legend()

kl77_ip25_pearson_coeff, kl77_ip25_pearson_p = pearsonr(kl77_ip25_centered, kl77_ngrip_centered)
kl77_ip25_spearman_coeff, kl77_ip25_spearman_p = spearmanr(kl77_ip25_centered, kl77_ngrip_centered)

# KL12 IP25 - plot centered IP25 values
kl12_ip25 = pd.read_csv(workdir+"input/interpol_IP25_SSTs_stack_average_praetorius_NAreplaced.txt", sep="\t")
kl12_ip25_ngrip = pd.merge_asof(kl12_ip25, ngrip_full, on='Age', direction='nearest')  # get nearest ngrip values to the timestamps of the ip25 values
kl12_ip25_ngrip.columns= ["SST", "IP25", "Age", "NGRIP"]

kl12_ip25 = kl12_ip25_ngrip["IP25"].values
kl12_ip25_centered = (kl12_ip25 - np.mean(kl12_ip25)) / np.std(kl12_ip25, ddof=0)
kl12_ngrip = kl12_ip25_ngrip["NGRIP"].values
kl12_ngrip_centered = (kl12_ngrip - np.mean(kl12_ngrip)) / np.std(kl12_ngrip, ddof=0)

plt.plot(kl12_ip25_ngrip["Age"]*1000, kl12_ip25_centered, label="KL12 IP$_{25}$", color="#17becf")
plt.plot(kl12_ip25_ngrip["Age"]*1000, kl12_ngrip_centered, label="NGRIP $\delta^{18}O$", color="#d62728")
plt.legend()
plt.title("Comparison of IP25 Records from Cores KL12, KL77 and NGRIP δ¹⁸O")
plt.xlabel("kyrs BP")
plt.ylabel("Standardized Proxy Values")
plt.savefig(workdir+"plots/ngrip_ip25.png")
plt.show()

# Correlate (spearman, pearson) IP25 from KL77 and KL12 with NGRIP
kl12_ip25_pearson_coeff, kl12_ip25_pearson_p = pearsonr(kl12_ip25_centered, kl12_ngrip_centered)
kl12_ip25_spearman_coeff, kl12_ip25_spearman_p = spearmanr(kl12_ip25_centered, kl12_ngrip_centered)

print("KL12 - IP25 | pearson coeff: ", round(kl12_ip25_pearson_coeff,2), " p: ", round(kl12_ip25_pearson_p,6))
print("KL12 - IP25 | spearman coeff: ", round(kl12_ip25_spearman_coeff,2), " p: ", round(kl12_ip25_spearman_p,6))

print("KL77 - IP25 | pearson coeff: ", round(kl77_ip25_pearson_coeff,2), " p: ", round(kl77_ip25_pearson_p,6))
print("KL77 - IP25 | spearman coeff: ", round(kl77_ip25_spearman_coeff,2), " p: ", round(kl77_ip25_spearman_p,6))