# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 11:26:28 2024

@author: Luisa
"""
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from scipy import stats
import re
sns.set(font="Arial")
plt.rcParams['pdf.fonttype']=42

df = pd.read_csv("20240506_pSTY_Site_annotation_Input.csv",sep=',')


df=df[['SiteID', 'com_dist', 'min_dist', 'cog_dist', 'res', 'DomainID','DomainID_seq', 'seq_wIndex', 'DomainNames', 'DomainType']]

# Calculate mean, range, standard deviation, and median for com_dist and min_dist for each SiteID and DomainID
result_df = df.groupby(['SiteID','DomainID','res','DomainType','DomainNames','DomainID_seq', 'seq_wIndex']).agg({
    'com_dist': ['mean', lambda x: x.max() - x.min(), 'std', 'median'],
    'cog_dist': ['mean', lambda x: x.max() - x.min(), 'std', 'median'],
    'min_dist': ['mean', lambda x: x.max() - x.min(), 'std', 'median']
}).reset_index()

# Rename columns for clarity
result_df.columns = ['SiteID','DomainID','res','DomainType','DomainNames','DomainID_seq', 'seq_wIndex', 'com_dist_mean', 'com_dist_range', 'com_dist_std', 'com_dist_median','cog_dist_mean', 'cog_dist_range', 'cog_dist_std', 'cog_dist_median',
                     'min_dist_mean', 'min_dist_range', 'min_dist_std', 'min_dist_median']



unique_Domains = result_df[['SiteID','DomainType','seq_wIndex']].drop_duplicates()
unique_Sites = result_df[['SiteID']].drop_duplicates()
unique_Sites[['UniprotID', 'Site']]=unique_Sites['SiteID'].str.split('_',expand=True)
unique_Enzymes=unique_Sites[['UniprotID']].drop_duplicates()



myPalette=sns.color_palette(['#E1CDF7','#80B1D3','#648FFF','#785EF0','#FD076D','#FE6100', '#8BE870','#1DE8C9'])
myPalette1=sns.color_palette(['#648FFF','#FE6100','#8BE870','#80B1D3','#785EF0','#FD076D','#1DE8C9','#E1CDF7'])
myPalette2={'S':'#1A85FF', 'Y':'#FFE56A', 'T':'#A5193E'}

sns.set_theme(style='white')
sns.histplot(data=result_df, x='com_dist_mean', bins=50, hue='res', multiple='stack', palette=myPalette2, kde=False)

# Calculate the 50th percentile
percentile_50 = result_df['com_dist_mean'].quantile(0.50)

# Add a vertical line at 75% of the data
plt.axvline(x=percentile_50, color='red', linestyle='dashed', linewidth=2, label='50th Percentile')
plt.tick_params(axis='both', which='both', left=True, bottom=True)
# Add text annotation for the red line
plt.text(percentile_50, plt.ylim()[1], f'  {percentile_50:.2f}', color='red', va='bottom', ha='left')

# Set labels and title
plt.xlabel('pSTY distance from domain COM (A)')
plt.ylabel('pSTY Count')

# Set x-axis limits to start from 0
plt.xlim(left=0)
plt.savefig('20240506_DomainCOM_histogram_50_res.pdf')


plot_df = result_df[~result_df['DomainType'].str.contains('Dimerization')]
plot_df = plot_df[~plot_df['DomainType'].str.contains('Other')]


# Get unique values from the column used for coloring and sort them
unique_categories = sorted(plot_df['DomainType'].unique())

myPalette=sns.color_palette(['#E1CDF7','#80B1D3','#648FFF','#785EF0','#FD076D','#FE6100', '#8BE870','#1DE8C9'])
# Define colors for each unique category
myPalette=sns.color_palette(['#FE6100','#785EF0','#FD076D','#648FFF','#E1CDF7','#8BE870'])
snspalette=sns.color_palette(palette='Set3')
colors = sns.color_palette(myPalette, len(unique_categories))

# Create a palette dictionary with unique categories as keys and assigned colors
palette_dict = {category: color for category, color in zip(unique_categories, colors)}

# Create a stacked histogram using Seaborn
sns.histplot(data=plot_df, x='res', hue='DomainType', multiple='stack', fill=True, hue_order=unique_categories,palette=palette_dict, binwidth=1, discrete=True, shrink=0.8, stat="count", legend=True,common_norm=False)
plt.xlabel('Amino Acid')
plt.ylabel('Frequency')

# Rotate x-axis ticks
# plt.xticks(rotation=90)
plt.yticks()
plt.tick_params(axis='both', which='both', left=True, bottom=True)

plt.savefig('20240506_Domains by res _histogram.pdf')
# plt.show()



















# Create a new DataFrame with rows containing substrings from the list in the DomainID column
filtered_df = result_df[result_df['com_dist_mean']<=23]
filtered_df = filtered_df[~filtered_df['DomainTypes'].str.contains('Dimerization')]
filtered_df = filtered_df[~filtered_df['DomainTypes'].str.contains('Other')]


Nucleotides=filtered_df[filtered_df['DomainTypes'].str.contains('Nucleotides')]
cofactors=filtered_df[filtered_df['DomainTypes'].str.contains('Cofactor')]
metals=filtered_df[filtered_df['DomainTypes'].str.contains('Metal')]
Actives=filtered_df[filtered_df['DomainTypes'].str.contains('Active')]
Substrates=filtered_df[filtered_df['DomainTypes'].str.contains('Substrate')]
Dimers=filtered_df[filtered_df['DomainTypes'].str.contains('Dimerization')]
catalytic=filtered_df[filtered_df['DomainTypes'].str.contains('Substrate') | filtered_df['DomainTypes'].str.contains('Active')]


result_df.to_csv('20240506_pSTY_domain_annotations_means.csv')
filtered_df.to_csv('20240506_pSTY_domain_annotations_23A.csv')




# Function to concatenate DomainTypes
def concatenate_domain_types(group):
    domain_types = ';'.join(sorted(set(group)))
    return domain_types

# Group by 'SiteID', 'Uniprot_seq_wIndex', 'res', and aggregate DomainTypes
new_alldf = alldf.groupby(['SiteID', 'Uniprot_seq_wIndex', 'res','com_dist_mean'])['DomainTypes'].agg(concatenate_domain_types).reset_index()

new_df23A=df_23A.groupby(['SiteID', 'Uniprot_seq_wIndex', 'res','com_dist_mean'])['DomainTypes'].agg(concatenate_domain_types).reset_index()


new_alldf.to_csv('20240502_pSTY_domain_annotations_all_noDup2.csv')
new_df23A.to_csv('20240502_pSTY_domain_annotations_23A_noDup2.csv')




sns.histplot(data=new_alldf, x='com_dist_mean', bins=50, hue='res', multiple='stack', palette=myPalette2, kde=False)

# Calculate the 50th percentile
percentile_50 = new_alldf['com_dist_mean'].quantile(0.50)

# Add a vertical line at 75% of the data
plt.axvline(x=percentile_50, color='red', linestyle='dashed', linewidth=2, label='50th Percentile')
plt.tick_params(axis='both', which='both', left=True, bottom=True)
# Add text annotation for the red line
plt.text(percentile_50, plt.ylim()[1], f'  {percentile_50:.2f}', color='red', va='bottom', ha='left')

# Set labels and title
plt.xlabel('pSTY distance from domain COM (A)')
plt.ylabel('pSTY Count')

# Set x-axis limits to start from 0
plt.xlim(left=0)






















custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", palette=myPalette)

sns.reset_orig()
sns.set_theme(style='white')


filtered_df['DomainTypes'] = filtered_df['DomainTypes'].str.lower()

# Get unique values from the column used for coloring and sort them
unique_categories = sorted(filtered_df['DomainTypes'].unique())

myPalette=sns.color_palette(['#E1CDF7','#80B1D3','#648FFF','#785EF0','#FD076D','#FE6100', '#8BE870','#1DE8C9'])
# Define colors for each unique category
myPalette=sns.color_palette(['#FE6100','#785EF0','#FD076D','#648FFF','#E1CDF7','#8BE870'])
colors = sns.color_palette(myPalette, len(unique_categories))

# Create a palette dictionary with unique categories as keys and assigned colors
palette_dict = {category: color for category, color in zip(unique_categories, colors)}

# Plot histogram for each category in Annotations
sns.histplot(data=filtered_df, x='com_dist_mean', bins=35, hue='DomainTypes', multiple='stack', fill=True, hue_order=unique_categories,palette=palette_dict,kde=False)
# plt.title('Dimerization Domain')
plt.xlabel('Distance from domain COM (Å)')
plt.ylabel('Frequency of STY')
plt.yticks()
plt.xlim(left=0)
plt.xlim(right=23)
plt.tick_params(axis='both', which='both', left=True, bottom=True)

plt.savefig('202405_STY by Domain COM_dist_histogram_23A.pdf')







# Get unique values from the column used for coloring and sort them
unique_categories = sorted(new_alldf['DomainTypes'].unique())

myPalette=sns.color_palette(['#E1CDF7','#80B1D3','#648FFF','#785EF0','#FD076D','#FE6100', '#8BE870','#1DE8C9'])
# Define colors for each unique category
myPalette=sns.color_palette(['#FE6100','#785EF0','#FD076D','#648FFF','#E1CDF7','#8BE870'])
snspalette=sns.color_palette(palette='Set3')
colors = sns.color_palette(snspalette, len(unique_categories))

# Create a palette dictionary with unique categories as keys and assigned colors
palette_dict = {category: color for category, color in zip(unique_categories, colors)}

# Create a stacked histogram using Seaborn
sns.histplot(data=new_alldf, x='res', hue='DomainTypes', multiple='stack', fill=True, hue_order=unique_categories,palette=palette_dict, binwidth=1, discrete=True, shrink=0.8, stat="count", legend=False,common_norm=False)
plt.xlabel('Amino Acid')
plt.ylabel('Frequency')

# Rotate x-axis ticks
# plt.xticks(rotation=90)
plt.yticks()
plt.tick_params(axis='both', which='both', left=True, bottom=True)

plt.savefig('Domains by res _histogram.pdf')
# plt.show()

unique_names = sorted(Actives['DomainNames'].unique())
top_values = Actives['DomainNames'].value_counts().sort_values(ascending=False).index[:5]

value_counts_df = Actives['DomainNames'].value_counts().reset_index()
value_counts_df.columns = ['DomainName', 'Frequency']
value_counts_df['Percentage'] = (value_counts_df['Frequency'] / value_counts_df['Frequency'].sum()) * 100


# Plotting KDE plot
sns.histplot(data=catalytic[catalytic['DomainNames'].map(catalytic['DomainNames'].value_counts()) > 50], x='com_dist_mean', hue='DomainNames', multiple='stack', palette=myPalette2, binwidth=1, discrete=True, shrink=0.8, stat="count", legend=True,common_norm=False)
# plt.xlim(left=0)
# plt.xlim(right=5)
# Set plot title and labels
plt.title('catalytic - com_dist_mean')
plt.xlabel('com_dist_mean')
plt.ylabel('pSTY Density')


testFilters=metals[metals['DomainNames'].map(metals['DomainNames'].value_counts()) > 20]
testFilters.loc[:, 'DomainNames'] = testFilters['DomainNames'].astype(str).apply(lambda x: x.lstrip('_'))


sns.histplot(data=testFilters, x='com_dist_mean', hue='DomainNames', multiple='stack')
# Get unique DomainNames and corresponding colors
unique_domains = testFilters['DomainID'].unique()
colors = sns.color_palette()[:len(unique_domains)]
# Create a dictionary mapping DomainNames to colors
legend_dict = {domain: color for domain, color in zip(unique_domains, colors)}

# Plot histplot with customized legend
sns.histplot(data=testFilters, x='com_dist_mean', hue='res', multiple='stack')




# Set a color palette
sns.set_palette('colorblind')

# Plotting KDE plot
sns.histplot(data=testFilters, x='com_dist_mean', multiple='fill',hue='DomainNames',common_norm=True)
# plt.xlim(left=0)
plt.xlim(right=25)
# Set plot title and labels
plt.title('KDE Plot of Mean com_dist_mean')
plt.xlabel('com_dist_mean')
plt.ylabel('pSTY Density')

# Add legend with hue_order
plt.legend(labels=filtered_df['res'].unique(), loc='upper right')

# Show the plot
plt.show()








# Dictionary of sidechain centers
sidechain_centers = {
    'ALA': ['CB'],
    'ARG': ['NE', 'NH1', 'NH2'],
    'ASN': ['OD1', 'ND2'],
    'ASP': ['OD1', 'OD2'],
    'CYS': ['SG'],
    'GLN': ['OE1', 'NE2'],
    'GLU': ['OE1', 'OE2'],
    'GLY': ['CA'],
    'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2'],
    'ILE': ['CD1'],
    'LEU': ['CD1', 'CD2'],
    'LYS': ['NZ'],
    'MET': ['SD'],
    'MSE': ['SE'],
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'PRO': ['CB', 'CG', 'CD'],
    'SER': ['OG'],
    'THR': ['OG1', 'CG2'],
    'TRP': ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3'],
    'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
    'VAL': ['CG1', 'CG2']
}

# Function to determine if atom is a sidechain center
def is_sidechain_center(residue, atom):
    if pd.isna(residue) or pd.isna(atom):
        return 0
    centers = sidechain_centers.get(residue, [])
    print(f"Residue: {residue}, Atom: {atom}, Centers: {centers}")  # Debug statement
    match = 1 if atom in centers else 0
    print(f"Match: {match}")  # Debug statement
    return match


# Function to determine if atom is a sidechain center
def keyLigands(residue):
    if residue in keyLigands:
        if residue in k and atom in v:
            print (residue, atom)
            return 1
        else:
            return 0



interface_dfr = pd.read_csv("PISA_r1r2_phos.csv",sep=',')
interface_dfrsub = interface_dfr[['uID','res1','atname1','res2','atname2', 'interface_pval','css','dist']]

# Ensure the columns are strings
interface_dfrsub.loc[:,'res1'] = interface_dfrsub['res1'].astype('string')
interface_dfrsub.loc[:,'atname1'] = interface_dfrsub['atname1'].astype('string')
interface_dfrsub.loc[:,'res2'] = interface_dfrsub['res2'].astype('string')
interface_dfrsub.loc[:,'atname2'] = interface_dfrsub['atname2'].astype('string')

# Strip whitespace from atom names
interface_dfrsub.loc[:,'atname1'] = interface_dfrsub['atname1'].str.strip()
interface_dfrsub.loc[:,'atname2'] = interface_dfrsub['atname2'].str.strip()

# Apply the function to create new columns
interface_dfrsub.loc[:,'res1_sidechain'] = interface_dfrsub.apply(lambda row: is_sidechain_center(row['res1'], row['atname1']), axis=1)
interface_dfrsub.loc[:,'res2_sidechain'] = interface_dfrsub.apply(lambda row: is_sidechain_center(row['res2'], row['atname2']), axis=1)



RaveInterface_df = interface_dfrsub.groupby(['uID','res1','res1_sidechain','res2','res2_sidechain']).agg({
    'interface_pval': ['mean', 'std', 'median','max','min', 'count'],
    'css': ['mean', 'std', 'median','max'],
    'dist': ['mean', 'std', 'median','max']
}).reset_index()

RaveInterface_df.columns = RaveInterface_df.columns.map('_'.join)
RfInterface_df=RaveInterface_df[RaveInterface_df['css_mean']>0.3]

# List of canonical amino acids
cAA = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

keyLigands=['2PG','3PG','F', 'K', 'U', 'AA', 'CA', 'CL', 'DA', 'DC', 'DG', 'DT', 'MG', 'MN', 'NA', 'SR', 'VX', 'YB', 'ZN', 'ABY', 'ACE', 
         'ACN', 'ACO', 'ACT', 'ACY', 'ADN', 'ADV', 'AGP', 'AHE', 'AKG', 'ALA', 'AMP', 'AMZ', 'AND', 'AOM', 'AOX', 'AOY', 'APO','APR', 'AQ1', 
         'ARE', 'ARG', 'ARH', 'ARS', 'ASD', 'ASN', 'ASO', 'ASP', 'AVL', 'AVQ', 'AVR', 'AVT', 'AYR', 'AZI', 'AZN', 'B12', 'B9D', 'BAU', 'BCO', 
         'BDP', 'BDT', 'BEF', 'BGC', 'BLA', 'BMP', 'BU2', 'BWS', 'C5P', 'CAA', 'CIT', 'CME', 'COA', 'COO', 'COS', 'CSO', 'CYS', 'DG2', 'DND', 
         'DTY', 'DUD', 'DUP', 'DUR', 'ENO', 'EOH','F6P', 'FAD', 'FBP', 'FLC', 'FMN', 'FUM', 'G1P', 'G6D', 'G6P', 'G6Q', 'GAI', 'GAR', 'GBP', 'GDN', 
         'GDP', 'GDS', 'GDU', 'GLC', 'GLN', 'GLU', 'GLV', 'GLY', 'GSH', 'GSN', 'GTX', 'GVX', 'H4B', 'HEM', 'HIS', 'HMG', 'HOH', 'HOX', 'ILE', 
         'IMD', 'IMP', 'IOD', 'IPR', 'ISN', 'LEU', 'LMR', 'LPA', 'LYS', 'M3L', 'MET', 'MLC', 'MLI', 'MLT', 'MLY', 'MLZ', 'MSE', 'NAD', 'NAG', 
         'NAI', 'NAP', 'NCA', 'NCN', 'NDP', 'NEP', 'NH2', 'NIO', 'NLE', 'NMN', 'NO3', 'NUP', 'O6A', 'OAA', 'OAS', 'OHP', 'OMP', 'ORO', 'OXL', 
         'OXM', 'PAH', 'PEP', 'PG2', 'PG4', 'PGA', 'PHE', 'PLP', 'PO3', 'PO4', 'POP', 'PPS', 'PPV', 'PRO', 'PRP', 'PUT', 'PYO', 'PYR', 'RBF', 
         'RBL', 'SAH', 'SAM', 'SCY', 'SEP', 'SER', 'SIN', 'SO4', 'STR', 'TAR', 'TES', 'THR', 'TPP', 'TRP', 'TTN', 'TUY', 'TVQ', 'TYR', 'TZ8', 
         'TZD', 'U1P', 'U5P', 'U91', 'UCY', 'UD1', 'UD2', 'UEP', 'UFP', 'UFT', 'UMP', 'UNU', 'UNX', 'UP6', 'UPG', 'URC', 'URE', 'UZT', 'VAL', 
         'VDN', 'VFV', 'VGN', 'VJJ', 'VU7', 'VWW', 'W7A', 'W8X', 'WDS', 'WDT', 'WDU', 'WDX', 'WW2', 'WWF', 'WY1', 'X74', 'X7J', 'X7P', 'X8P', 
         'XMP', 'XSP', 'XV1', 'XX0', 'Y0X', 'Y9B', 'YBY', 'YF7', 'YOJ', 'YQP', 'YY1', 'YY2', 'YZ9', 'ZEC', 'ZES', 'ZIN', 'ZOM', 'ZON', 'ZST', 'ZXU', 'ZXY']


# Filter rows where 'res2' is a canonical amino acid
aaInterface_df = RaveInterface_df[RaveInterface_df['res2_'].isin(cAA)]

notaaInterface_df = RaveInterface_df[~RaveInterface_df['res2_'].isin(cAA)]


keyR2Inteface_df=notaaInterface_df[notaaInterface_df['res2_'].isin(keyLigands)]


sns.reset_orig()
sns.set_theme(style='ticks')
# Plot histogram for each category in Annotations
sns.histplot(data=keyR2Inteface_df, x='css_mean', bins=20, hue='res1_', hue_order=sorted(keyR2Inteface_df['res1_'].unique()) ,multiple='stack', palette=myPalette3,kde=False)
plt.title('Histogram of interface Score')
plt.xlabel('Probability of role in dimer interface (CSS score)')
plt.ylabel('Frequency')
plt.xlim(left=0)
plt.xlim(right=1)
plt.yticks()
plt.tick_params(axis='both', which='both', left=True, bottom=True)
plt.savefig('20240521_pSTY_interface_css03.pdf')


RaveInterface_df.to_csv("20240726_Interface_CSS_averaged.csv")
RfInterface_df.to_csv("20240726_Interface_CSS_over03.csv")

aaInterface_df.to_csv("20240726_Interface_aaInteractors.csv")
keyR2Inteface_df.to_csv("20240726_Interface_notaaInteractors.csv")


interface_df = pd.read_csv("filtered_phos_residues.csv",sep=',')
interface_dfsub = interface_df[['uID','name','pvalue','solv_en_over_bsa','interface_pval', 'interface_solv_en','stab_en','css', 'interface_type']]

# Define a function to apply multiple hypothesis correction
def correct_pvalues(df):
    corrected_pvals = multipletests(df['pvalue'], method='fdr_bh')[1]
    df['corrected_pvalue'] = corrected_pvals
    return df

# Group by 'uID' and apply the correction
corrected_df = interface_dfsub.groupby('uID').apply(correct_pvalues).reset_index(drop=True)


aveInterface_df = interface_dfsub.groupby(['uID','name','interface_type']).agg({
    'pvalue': ['mean', 'std', 'median','max','min'],
    'solv_en_over_bsa': ['mean', 'std', 'median','max'],
    'interface_pval': ['mean', 'std', 'median','max'],
    'interface_solv_en': ['mean', 'std', 'median','max'],
    'stab_en': ['mean', 'std', 'median','max'],
    'css': ['mean', 'std', 'median','max']
}).reset_index()

aveInterface_df.columns = aveInterface_df.columns.map('_'.join)
fInterface_df=aveInterface_df[aveInterface_df['css_max']>0.3]

aveInterface_df.to_csv("Interface_CSS_averaged.csv")
fInterface_df.to_csv("Interface_CSS_over03.csv")



sns.reset_orig()
sns.set_theme(style='ticks')
# Plot histogram for each category in Annotations
sns.histplot(data=aveInterface_df, x='css_max', bins=20, hue='name_', hue_order=sorted(aveInterface_df['name_'].unique()) ,multiple='stack', palette=myPalette3,kde=False)
plt.title('Histogram of interface Score')
plt.xlabel('Probability of role in dimer interface (CSS score)')
plt.ylabel('Frequency')
plt.xlim(left=0)
plt.xlim(right=1)
plt.yticks()
plt.tick_params(axis='both', which='both', left=True, bottom=True)
plt.savefig('20240521_pSTY_interface_css03.pdf')


sns.regplot(x='pvalue_mean', y='css_mean', data=aveInterface_df, line_kws=myPalette3)
# plt.title(f'Pearson Correlation: {correlation:.2f}')
plt.xlabel('pvalue_mean')
plt.ylabel('css_mean')
plt.show()

sns.pairplot(data=interface_dfsub[['pvalue','css','name']], hue='name',palette=myPalette3)
# plt.title(f'Pearson Correlation: {correlation:.2f}')
plt.xlabel('pvalue_mean')
plt.ylabel('css_mean')
plt.show()



# Count the frequency of each amino acid
aaCounts = aaInterface_df['res2_'].value_counts()

# Sort res2 categories by frequency
sorted_res2 = aaCounts.index

# Sort the DataFrame by res2 frequency
aaInterface_df.loc[:, 'res2_'] = pd.Categorical(aaInterface_df['res2_'], categories=sorted_res2, ordered=True)
aaInterface_df = aaInterface_df.sort_values('res2_')

# Apply the sorted order to the 'res2_' column in the DataFrame
aaInterface_df['res2_'] = pd.Categorical(aaInterface_df['res2_'], categories=sorted_res2, ordered=True)



# Plotting
plt.bar(aaCounts.index, aaCounts.values, color='skyblue')
plt.xlabel('Amino Acid')
plt.ylabel('Frequency')
# Rotate x-axis ticks
plt.xticks(rotation=90)
plt.title('Frequency of Amino Acids')
plt.show()


myPalette2=sns.set_palette(sns.color_palette(['#6a3d9a','#1f78b4','#33a02c']))

# Create a stacked histogram using Seaborn
sns.histplot(data=aaInterface_df, x='res2_', hue='res1_', multiple='fill', palette=myPalette3, binwidth=1, discrete=True, shrink=0.8, stat="count", common_norm=False)
plt.xlabel('Amino Acid')
plt.ylabel('Density')

# Rotate x-axis ticks
plt.xticks(rotation=90)
plt.yticks()
plt.tick_params(axis='both', which='both', left=True, bottom=True)

sns.set(font="Arial")
plt.rcParams['pdf.fonttype']=42
sns.reset_orig()
sns.set_theme(style='ticks')
sns.set(font="Arial")
plt.rcParams['pdf.fonttype']=42


# Plotting KDE plot
sns.kdeplot(data=RaveInterface_df, x='dist_mean', fill=True, hue='res1_', hue_order=sorted(RaveInterface_df['res1_'].unique()), common_norm=True,palette=myPalette3)

# Set plot title and labels
plt.title('KDE Plot of Mean distance')
plt.xlabel('Average distance between interactions (Å)')
plt.ylabel('pSTY Density')


plt.savefig('pSTY_interface_meanDist_kde.pdf')



# Create a DataFrame with counts
aaInterface_df['res2_wSidechain']=aaInterface_df['res2_']+'_'+aaInterface_df['res2_sidechain_'].astype('string')
aaInterface_df['res1_wSidechain']=aaInterface_df['res1_']+'_'+aaInterface_df['res1_sidechain_'].astype('string')
counts_df = aaInterface_df.groupby(['res1_wSidechain','res2_wSidechain']).size().reset_index(name='count')


# Pivot the DataFrame to create a stacked bar plot
pivot_df = counts_df.pivot(index=['res2_wSidechain'], columns=['res1_wSidechain'], values='count')


# Fill NaN values with 0 (or another appropriate value)
pivot_df = pivot_df.fillna(0)

# Create the clustered heatmap using Seaborn's clustermap
g = sns.clustermap(pivot_df, annot=False, cmap = sns.light_palette("seagreen", 15, as_cmap=False),vmin=0, center=15,vmax=np.max(pivot_df),metric='euclidean', method='complete',cbar=True, row_cluster=True, col_cluster=True,yticklabels=True,xticklabels=True ,annot_kws={"color": "black"}, figsize=(5, 10))
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)

plt.savefig('pSTY_interface_notaaCount_heatmapT.pdf')


sns.reset_orig()
sns.set_theme(style='ticks')

# Plotting KDE plot
sns.kdeplot(data=aaInterface_df, x='dist_mean', fill=True, hue='res1_', hue_order=sorted(aaInterface_df['res1_'].unique()) ,common_norm=True,palette=myPalette3)

# Set plot title and labels
plt.title('KDE Plot of Mean distance')
plt.xlabel('Average distance between interactions (Å)')
plt.ylabel('pSTY Density')

# Add legend with hue_order
plt.legend(labels=aaInterface_df['res1_'].unique(), loc='upper right')

# Show the plot
plt.show()


plt.savefig('pSTY_interface_kde.pdf')




# Create a DataFrame with counts
keyR2Inteface_df['res2_wSidechain']=keyR2Inteface_df['res2_']+'_'+keyR2Inteface_df['res2_sidechain_'].astype('string')
keyR2Inteface_df['res1_wSidechain']=keyR2Inteface_df['res1_']+'_'+keyR2Inteface_df['res1_sidechain_'].astype('string')
counts_df = keyR2Inteface_df.groupby(['res1_wSidechain','res2_']).size().reset_index(name='count')


# Pivot the DataFrame to create a stacked bar plot
pivot_df = counts_df.pivot(index=['res2_'], columns=['res1_wSidechain'], values='count')


# Fill NaN values with 0 (or another appropriate value)
pivot_df = pivot_df.fillna(0)

# Create the clustered heatmap using Seaborn's clustermap
g = sns.clustermap(pivot_df, annot=False, cmap = sns.light_palette("seagreen", 15, as_cmap=False),vmax=np.max(pivot_df),metric='euclidean', method='complete',cbar=True, row_cluster=True, col_cluster=True,yticklabels=True,xticklabels=True ,annot_kws={"color": "black"}, figsize=(5, 10))
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)

plt.savefig('pSTY_interface_notaaCount_heatmapT.pdf')

# Create a DataFrame with counts
notaaInterface_df['res1_wSidechain']=keyR2Inteface_df['res1_']+'_'+keyR2Inteface_df['res1_sidechain_'].astype('string')
counts_df = notaaInterface_df.groupby(['res2_', 'res1_wSidechain']).size().reset_index(name='count')


# Pivot the DataFrame to create a stacked bar plot
pivot_df = counts_df.pivot(index='res2_', columns=['res1_wSidechain'], values='count')



# Fill NaN values with 0 (or another appropriate value)
pivot_df = pivot_df.fillna(0)

# Create the clustered heatmap using Seaborn's clustermap
g = sns.clustermap(pivot_df, annot=False, cmap = sns.light_palette("seagreen", 15, as_cmap=False),vmax=np.max(pivot_df),metric='euclidean', method='complete',cbar=True, row_cluster=True, col_cluster=True,yticklabels=True,xticklabels=True ,annot_kws={"color": "black"}, figsize=(5, 10))
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
# plt.show()

plt.savefig('pSTY_interface_aa_heatmap.pdf')

myPalette3={'SER':'#1A85FF', 'THR':'#A5193E', 'TYR':'#FFB000'}


# Set the order of columns based on the custom color palette
column_order = [ 'SER','THR','TYR']

# Reorder the columns in the DataFrame
pivot_df = pivot_df[column_order]

# Plot the stacked bar plot
plt.figure(figsize=(10, 6))
pivot_df.plot(kind='bar', stacked=True, color=myPalette3.values(), edgecolor='white', width=0.8)


plt.xlabel('Amino Acid')
plt.ylabel('Percentage of Total')
plt.title('Percentage of Amino Acid Distribution by res1')

# Rotate x-axis ticks
plt.xticks(rotation=90)


# plt.show()
plt.savefig('pSTY_interface_aa_dist.pdf')









# Create a DataFrame with counts
counts_doms = filtered_df.groupby(['res', 'Domains']).size().reset_index(name='count')

# Calculate the percentage over the entire dataset
counts_doms['percentage'] = counts_doms.groupby('Domains')['count'].transform(lambda x: x / x.sum() * 100)
counts_doms['percentage'] = counts_doms['count'] / counts_doms['count'].sum() * 100

# Pivot the DataFrame to create a stacked bar plot
pivot_doms = counts_doms.pivot(index='Domains', columns='res', values='percentage')


myPalette4={'S':'#1A85FF', 'T':'#A5193E', 'Y':'#FFB000'}
myPalette5=sns.color_palette(['#FFB000','#BEBADA','#80B1D3','#648FFF','#785EF0','#FD076D','#FE6100'])

# Set the order of columns based on the custom color palette
index_order = [ 'S','Y','T']

# Reorder the index of the DataFrame
pivot_doms = pivot_doms.reindex(index_order)


pivot_doms.plot(kind='bar', stacked=True, color=myPalette5[::-1], edgecolor='white')

plt.ylabel('Percentage of Total')
plt.title('Percentage of STY by domain')


plt.legend(loc='upper left', bbox_to_anchor=(1, 1))


# plt.show()
plt.savefig('Percentage of STY by domain.pdf')














# Set a color palette
sns.set_palette('colorblind')

# Plotting KDE plot
sns.kdeplot(data=fInterface_df, x='dist_mean_', fill=True, hue='res1_', common_norm=True)

# Set plot title and labels
plt.title('KDE Plot of Mean css_mean')
plt.xlabel('CSS Score')
plt.ylabel('pSTY Density')

# Add legend with hue_order
plt.legend(labels=fInterface_df['res1_'].unique(), loc='upper right')

# Show the plot
plt.show()


plt.savefig('pSTY_interface_kde.pdf')










# Count the occurrences of each domain
domain_counts = result_df['Domains'].value_counts()

# Create a pie chart
plt.figure(figsize=(8, 8))
plt.pie(domain_counts, labels=domain_counts.index, autopct='%1.1f%%', startangle=140, colors=sns.color_palette('pastel'))
plt.title('Distribution of Domains')
plt.show()





# Plot histogram for each category in Annotations
sns.histplot(data=filtered_df, x='com_dist_mean', bins=100, hue='Domains', multiple='stack', palette="bright")
plt.title('Distribution of pSTY distance from COM of Domains')
plt.xlabel('Mean distance to COM (A)')
plt.ylabel('pSTY Frequency')
plt.show()
plt.clf()



# Set up the figure
plt.figure(figsize=(12, 8))

# Create a KDE plot for each category in the 'Domains' column
sns.histplot(
    data=filtered_df,
    x='com_dist_mean',
    hue='Domains',
    multiple='fill',
    fill=True,
    common_norm=False,
    palette='colorblind',
    alpha=0.5,
    linewidth=0,
    kde=True
)



# Add labels and title
plt.xlabel('com_dist_mean')
plt.ylabel('Density')
plt.title('KDE Plot of com_dist_mean for Different Domains')

# Show the plot
plt.show()





# Calculate Pearson correlation
correlation = filtered_df['com_dist_mean'].corr(filtered_df['min_dist_mean'])

# Plot the scatter plot with regression line
sns.regplot(x='com_dist_mean', y='min_dist_mean', data=filtered_df)
plt.title(f'Pearson Correlation: {correlation:.2f}')
plt.xlabel('com_dist_mean')
plt.ylabel('min_dist_mean')
plt.show()




# Assuming result_df is the DataFrame containing the means
# Replace result_df with the actual variable name if different

# Plotting histogram
plt.hist(result_df['com_dist_mean'], bins=50, edgecolor='black')
plt.title('Histogram of Mean COM distance from domain of interest')
plt.xlabel('Distance in Angstromes')
plt.ylabel('Frequency')
plt.show()




# Assuming result_df is the DataFrame containing the means
# Replace result_df with the actual variable name if different


# Sorting DataFrame by com_dist_mean
result_df = result_df.sort_values('com_dist_mean')

# Creating a waterfall plot
plt.bar(result_df['SiteID'] + ' ' + result_df['DomainID'], result_df['com_dist_mean'], color='blue')
plt.title('Waterfall Plot of Mean com_dist_mean')
plt.xlabel('SiteID - DomainID')
plt.ylabel('Mean com_dist_mean')
plt.xticks([])
plt.show()







df_23A = pd.read_csv("STY_ActiveSite_23A.csv",sep=',')
df_12A = pd.read_csv("20240518_STY_count_12A.csv",sep=',')

prot_12A=df_12A.drop(['UniprotDom','sequence','domainName','num_SER','num_THR', 'num_TYR'], axis=1)
prot_12A['UniprotID']=prot_12A['UniprotID'].str.split('_', expand=True)[0]



df=prot_12A.copy()

import ast

def ensure_list(column):
    return column.apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

df['SER_index'] = ensure_list(df['SER_index'])
df['THR_index'] = ensure_list(df['THR_index'])
df['TYR_index'] = ensure_list(df['TYR_index'])


def combine_and_deduplicate_lists(series):
    combined_list = [item for sublist in series for item in sublist]
    deduplicated_list = list(set(combined_list))  # Removes duplicates
    deduplicated_list.sort(key=combined_list.index)  # Preserve the original order of the first occurrence
    return deduplicated_list


# Using groupby and agg with the custom function
combined_df = df.groupby(['UniprotID','UniprotDom','domainName','sequence','UniprotSeq'], as_index=False).agg({
    'SER_index': combine_and_deduplicate_lists,
    'THR_index': combine_and_deduplicate_lists,
    'TYR_index': combine_and_deduplicate_lists
})


# Using groupby and agg with the custom function
combined_dfprot = df.groupby(['UniprotID','UniprotSeq'], as_index=False).agg({
    'SER_index': combine_and_deduplicate_lists,
    'THR_index': combine_and_deduplicate_lists,
    'TYR_index': combine_and_deduplicate_lists
})

# Using groupby and agg with the custom function
combined_dom = df.groupby('UniprotID', as_index=False).agg({
    'SER_index': combine_and_deduplicate_lists,
    'THR_index': combine_and_deduplicate_lists,
    'TYR_index': combine_and_deduplicate_lists
})




df2=combined_dfprot.copy()

def verify_and_correct_indices(row):
    uniprot_seq = row['UniprotSeq']
    
    def correct_index(residue_list, residue_char):
        corrected_indices = set()  # Use set to avoid duplicates
        for chain, residue, idx in residue_list:
            idx = int(idx) - 1  # Convert to 0-based index
            if 0 <= idx < len(uniprot_seq) and uniprot_seq[idx] == residue_char:
                corrected_indices.add((chain, residue, str(idx + 1)))
            elif 0 <= (idx + 1) < len(uniprot_seq) and uniprot_seq[idx + 1] == residue_char:
                corrected_indices.add((chain, residue, str(idx + 2)))  # Correct off-by-one error
            # If neither the original nor the +1 index is correct, do not add to corrected_indices
        return corrected_indices  # Keep as set to maintain unique elements

    row['SER_index'] = correct_index(row['SER_index'], 'S')
    row['THR_index'] = correct_index(row['THR_index'], 'T')
    row['TYR_index'] = correct_index(row['TYR_index'], 'Y')
    
    return row

# Apply index correction
df2 = df2.apply(verify_and_correct_indices, axis=1)



# Remove chain ID and duplicates
def remove_chain_and_duplicates(residue_set):
    return list({(residue, idx) for _, residue, idx in residue_set})

df2['SER_index'] = df2['SER_index'].apply(remove_chain_and_duplicates)
df2['THR_index'] = df2['THR_index'].apply(remove_chain_and_duplicates)
df2['TYR_index'] = df2['TYR_index'].apply(remove_chain_and_duplicates)


df2['num_SER'] = df2['SER_index'].apply(len)
df2['num_THR'] = df2['THR_index'].apply(len)
df2['num_TYR'] = df2['TYR_index'].apply(len)
df2['dom_STY'] = df2['num_SER'] + df2['num_THR'] + df2['num_TYR']

# Add columns for total counts of S, T, Y, and STY
df2['total_S'] = df2['UniprotSeq'].str.count('S')
df2['total_T'] = df2['UniprotSeq'].str.count('T')
df2['total_Y'] = df2['UniprotSeq'].str.count('Y')
df2['total_STY'] = df2['total_S'] + df2['total_T'] + df2['total_Y']

df2.to_csv('STY_count_12A_perProt.csv')


# Define a custom aggregation function to calculate the maximum value for each column
def max_values_per_group(group):
    max_values = group[['num_SER', 'num_THR', 'num_TYR']].max()
    return max_values


# Group by 'DomainID' and 'Seq_wIndex', then apply the custom aggregation function
max_23A = df_23A.groupby(['UniprotID', 'domainName', 'UniprotDom','UniprotSeq']).apply(max_values_per_group)
max_12A = df_12A.groupby(['UniprotID', 'DomainName']).apply(max_values_per_group)

# Preprocess data to merge rows with different 'DomainID' but same 'Seq_wIndex'
merged_df = df_23A.groupby(['UniprotID']).apply(lambda x: x.groupby('DomainName').max()).reset_index(drop=True)

# Group by 'Seq_wIndex' again and apply the custom aggregation function
max_values2 = merged_df.groupby(['seq_wIndex']).apply(max_values_per_group)

max_23A.to_csv('STY_Counts_23A_max.csv')
max_12A.to_csv('STY_Counts_12A_max.csv')

