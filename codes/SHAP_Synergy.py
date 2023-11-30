# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 12:12:31 2022

@author: PI
"""

import shap
from MLP import shap_values, PI_cutoff
import matplotlib.pyplot as plt
import itertools

# Selected features for plotting
selection_9 = ['qed', 'SMR_VSA4', 'RingCount', 'SlogP_VSA4', 'fr_Ndealkylation2', 'EState_VSA7', 'NumAromaticCarbocycles', 'fr_ether', 'PEOE_VSA1']

# Use itertools.combinations() to find all combinations of two elements
combinations = list(itertools.combinations(selection_9, 2))

# Counter for figure naming
i = 0

# Loop through all possible combinations
for combination in combinations:
    i += 1
    shap_list = list(combination)

    # Create dependence plot
    shap.dependence_plot(shap_list[0], shap_values, PI_cutoff.iloc[:, 1:][selection_9], interaction_index=shap_list[1], show=False)
    ax = plt.gca()

    # Customize plot appearance
    ax.spines['right'].set_visible(True)
    ax.spines['top'].set_visible(True)
    ax.spines['left'].set_linewidth(3)
    ax.spines['right'].set_linewidth(3)
    ax.spines['bottom'].set_linewidth(3)
    ax.spines['top'].set_linewidth(3)

    ax.set_xlabel(shap_list[0], fontname='Arial', fontsize=16)
    ax.set_ylabel(f'SHAP value for {shap_list[0]}', fontname='Arial', fontsize=16)

    ax.tick_params(axis='x', which='both', width=2)
    ax.tick_params(axis='y', which='both', width=2)

    x_ticklabels = ax.get_xticklabels()
    for label in x_ticklabels:
        label.set_fontproperties('Arial')
        label.set_fontsize(14)

    y_ticklabels = ax.get_yticklabels()
    for label in y_ticklabels:
        label.set_fontproperties('Arial')
        label.set_fontsize(14)

    # Save the figure with a unique name
    figure_name = f"{i}_{shap_list[0]}_{shap_list[1]}"
    plt.savefig(f'..\figures\{figure_name}.png', dpi=900)
