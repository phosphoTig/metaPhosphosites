# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 00:41:51 2023

@author: Luisa
"""
import pandas as pd
import numpy as np
from scipy.stats import hypergeom, chi2_contingency
from statsmodels.stats.multitest import multipletests

df=pd.read_csv('Hyperg_input_23A_Dom.csv', sep=',').set_index('DomainType')


# Define a function to perform the hypergeometric test for each row and calculate p-values
def hypergeom_test(row):
    k = row['k (phospho count w/in 23A)']  # number of successes in the sample
    n = row['s ( total count w/in 23A)']   # number of draws (sample size)
    M = row['N (total AA count)']          # total population size
    N = row['M (total phospho count)']     # number of successes in the population
    
    # Calculate the one-tailed p-value
    p_value_one_tailed = 1 - hypergeom.cdf(k - 1, M, N, n)
    
    return p_value_one_tailed

# Apply the function to each row of the DataFrame and store the results in new columns
df['p_value'] = df.apply(hypergeom_test, axis=1)


# Define the threshold for statistical significance
alpha = 0.05

# Correct for multiple hypothesis testing using the Benjamini-Hochberg procedure
reject, p_values_corrected, _, _ = multipletests(df['p_value'], alpha=alpha, method='fdr_bh')


df['FDR']=p_values_corrected
df['reject H0']=reject


df.to_csv("20240519_hyperGtset_23A.csv")


import matplotlib.pyplot as plt
import seaborn as sns
sns.set(font="Arial")
plt.rcParams['pdf.fonttype']=42
sns.set_theme(style='white')
# Plot a bar plot of p-values for each domain type
plt.figure(figsize=(6, 3))
df_sorted = df.sort_values('FDR')
sns.barplot(x=df_sorted.index, y=-np.log10(df_sorted['FDR']))
plt.xticks(rotation=90)
plt.yticks()
plt.tick_params(axis='both', which='both', left=True, bottom=True)
plt.title('Enrichment of sites < 23Å from Domain')
plt.xlabel('Domain Type')
plt.ylabel('-log10(FDR)')
plt.axhline(y=-np.log10(0.05), color='r', linestyle='--')  # Add a horizontal line at p=0.05
# plt.show()

plt.savefig('pSTY_hyperG_plot_23A.pdf')
