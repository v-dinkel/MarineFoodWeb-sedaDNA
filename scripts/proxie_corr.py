# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 10:15:47 2025

@author: vdinkel
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

# NGRIP from manuscript
ngrip = pd.read_table("C:/Users/vdinkel/Desktop/Manuscript/output/KL-77_ngrip.csv", sep=";", header=0)

# Full NGRIP from source
ngrip_full = pd.read_table("C:/Users/vdinkel/Desktop/Manuscript/input/ngrip.csv")
ngrip_full.columns = ["Age", "Temp"]

# KL77 IP25
kl77_ip25 = pd.read_csv("C:/Users/vdinkel/Desktop/Manuscript/submission_Nature_Communications/revision1/script/SST/IP25_77KL.csv", sep=";")
kl77_ip25_ngrip = pd.merge_asof(kl77_ip25, ngrip_full, on='Age', direction='nearest')
kl77_ip25_ngrip.columns = ["Age", "IP25", "NGRIP"]

kl77_ip25 = kl77_ip25_ngrip["IP25"].values
kl77_ip25_centered = (kl77_ip25 - np.mean(kl77_ip25)) / np.std(kl77_ip25, ddof=0)  # ddof=0 for population std
kl77_ngrip = kl77_ip25_ngrip["NGRIP"].values
kl77_ngrip_centered = (kl77_ngrip - np.mean(kl77_ngrip)) / np.std(kl77_ngrip, ddof=0)  # ddof=0 for population std

plt.figure(figsize=(8, 6), dpi=300)
plt.plot(kl77_ip25_ngrip["Age"]*1000, kl77_ip25_centered, label="KL77 IP$_{25}$", color="#1f77b4")
plt.legend()

kl77_ip25_pearson_coeff, kl77_ip25_pearson_p = pearsonr(kl77_ip25_centered, kl77_ngrip_centered)
kl77_ip25_spearman_coeff, kl77_ip25_spearman_p = spearmanr(kl77_ip25_centered, kl77_ngrip_centered)

# KL12 IP25
kl12_ip25 = pd.read_csv("C:/Users/vdinkel/Desktop/Manuscript/submission_Nature_Communications/revision1/script/SST/interpol_IP25_SSTs_stack_average_praetorius_NAreplaced.txt", sep="\t")
kl12_ip25_ngrip = pd.merge_asof(kl12_ip25, ngrip_full, on='Age', direction='nearest')
kl12_ip25_ngrip.columns= ["SST", "IP25", "Age", "NGRIP"]

kl12_ip25 = kl12_ip25_ngrip["IP25"].values
kl12_ip25_centered = (kl12_ip25 - np.mean(kl12_ip25)) / np.std(kl12_ip25, ddof=0)  # ddof=0 for population std
kl12_ngrip = kl12_ip25_ngrip["NGRIP"].values
kl12_ngrip_centered = (kl12_ngrip - np.mean(kl12_ngrip)) / np.std(kl12_ngrip, ddof=0)  # ddof=0 for population std

plt.plot(kl12_ip25_ngrip["Age"]*1000, kl12_ip25_centered, label="KL12 IP$_{25}$", color="#17becf")
plt.plot(kl12_ip25_ngrip["Age"]*1000, kl12_ngrip_centered, label="NGRIP $\delta^{18}O$", color="#d62728")
plt.legend()
plt.title("Comparison of IP25 Records from Cores KL12, KL77 and NGRIP δ¹⁸O")
plt.xlabel("kyrs BP")
plt.ylabel("Standardized Proxy Values")
plt.savefig("C:/Users/vdinkel/Desktop/Manuscript/submission_Nature_Communications/revision1/script/new_plots/ngrip_ip25.png")
plt.show()

kl12_ip25_pearson_coeff, kl12_ip25_pearson_p = pearsonr(kl12_ip25_centered, kl12_ngrip_centered)
kl12_ip25_spearman_coeff, kl12_ip25_spearman_p = spearmanr(kl12_ip25_centered, kl12_ngrip_centered)

print("KL12 - IP25 | pearson coeff: ", round(kl12_ip25_pearson_coeff,2), " p: ", round(kl12_ip25_pearson_p,6))
print("KL12 - IP25 | spearman coeff: ", round(kl12_ip25_spearman_coeff,2), " p: ", round(kl12_ip25_spearman_p,6))

print("KL77 - IP25 | pearson coeff: ", round(kl77_ip25_pearson_coeff,2), " p: ", round(kl77_ip25_pearson_p,6))
print("KL77 - IP25 | spearman coeff: ", round(kl77_ip25_spearman_coeff,2), " p: ", round(kl77_ip25_spearman_p,6))