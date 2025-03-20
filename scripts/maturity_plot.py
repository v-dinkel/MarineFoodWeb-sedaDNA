# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 16:01:43 2025

@author: vdinkel
"""

import matplotlib.pyplot as plt
import numpy as np

inc_col = "#333333" # orange ##D4A017
dec_col = "#CCCCCC" # #1F77B4

# Data
categories = ["Bottom-Up", "Top-Down", "Intraguild"]
glacial_values = [0.15, 0.15, 0.34]
interglacial_values = [0.08, 0.22, 0.38]

x = np.arange(len(categories))
width = 0.35  # Width of bars

fig, ax = plt.subplots(figsize=(7, 5))

# Plot bars for Glacial and Interglacial
ax.bar(x - width/2, glacial_values, width, label="IG2 Cluster 0", color=dec_col)
ax.bar(x + width/2, interglacial_values, width, label="IG2 Cluster 1", color=inc_col)

# Labels and legend
ax.set_ylabel("% of Relative Ascendency")
ax.set_xlabel("Interaction Type")
ax.set_title("Changes in IG1-Module Structure")
ax.set_xticks(x)
ax.set_xticklabels(categories)
ax.legend()

plt.savefig("C:/Users/vdinkel/Desktop/Manuscript/submission_Nature_Communications/revision1/script/new_plots/S7_bo_td_overview.png")
plt.show()