# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 10:54:45 2022

@author: Luisa
"""

import pandas as pd
import seaborn as sns
import plotly.graph_objects as go

import matplotlib.pyplot as plt
import plotly.figure_factory as ff

import missingno as msno

import plotly.io as pio
import plotly.express as px
from PIL import Image

pio.renderers.default = "svg"

sns.set(font="Arial")
plt.rcParams['pdf.fonttype']=42
sns.set_theme(style="ticks")

Res=pd.read_csv('2024_PDB_rvals.csv', sep=',').set_index('PDB_id')
Res=Res.dropna()

sns.distplot(Res['Resolution'], kde=True, vertical=True, rug=True, hist=False, kde_kws=dict(shade=True), rug_kws=dict(lw=2, color='orange'))
plt.savefig('Resolution PDB.pdf')




sns.distplot(Res['R_free - R_work'], kde=True, vertical=True, rug=True, hist=False, kde_kws=dict(shade=True), rug_kws=dict(lw=2, color='orange'))
plt.savefig("Density histogram R v Rfree.pdf")

sns.jointplot(x=Res['Resolution'], y=Res['R_free - R_work'], kind='reg',color='#80cdc1', xlim=(0,5),ylim=(-0.1,0.15))


sns.jointplot(x=Res['Resolution'], y=Res['R_free - R_work'], kind="kde", color='#80cdc1', xlim=(0,5),ylim=(-0.1,0.15), stat_func=None, marginal_kws={"color":"black", "lw":0.5}, joint_kws={"colors":"black","cmap":None, "linewidths":0.5}, shade=False, n_levels=10)
plt.savefig("R v Rfree Jointplot.pdf",bbox_inches='tight')
plt.clf()

Res=Res[Res['Resolution']<4]


fig = go.Figure()

fig.add_trace(go.Histogram2dContour(
    x=Res['Resolution'], y=Res['R_free - R_work'],
    xaxis = 'x', yaxis= 'y',
    colorscale = 'greens'
))


# Side dist
fig.add_trace(go.Histogram(
       y=Res['R_free - R_work'],
       xaxis = 'x2',
       marker = dict(color = 'rgba(0,0,0,1)'),
       ybins=dict(start=-0.05,end=max(Res['R_free - R_work'])),
       nbinsy=300
))

# Side dist
fig.add_trace(go.Histogram(
    x=Res['Resolution'],
    yaxis='y2',
    marker=dict(color='rgba(0,0,0,1)'),
    xbins=dict(start=min(Res['Resolution']), end=4),
    nbinsx=100
))


fig.update_layout(
    autosize=True,
    font=dict(family="Arial"),
    xaxis=dict(
        zeroline=True,
        showline=True,  # Set to False to remove the x-axis line
        domain=[0, 0.85],
        showgrid=False,
        tickmode='array',  # Set tick mode to 'array'
        tickvals=[1, 2, 3, 4],  # Set the tick values
        ticktext=['1', '2', '3', '4'],  # Set the tick labels
        ticks='outside',  # Set tick placement to 'outside'
        tickfont=dict(color='black'),
        range=[min(Res['Resolution']), 4],
        title='X Axis Title'  # Set x-axis title
    ),
    yaxis=dict(
        zeroline=False,
        showline=True,  # Set to False to remove the y-axis line
        domain=[0, 0.85],
        showgrid=False,
        tickmode='array',  # Set tick mode to 'array'
        tickvals=[-0.05, 0, 0.05, 0.1, 0.15],  # Set the tick values
        ticktext=[ '-0.05', '0', '0.05', '0.1', '0.15'],  # Set the tick labels
        ticks='outside',  # Set tick placement to 'outside'
        tickfont=dict(color='black'),
        range=[-0.05, max(Res['R_free - R_work'])],
        title='Y Axis Title'  # Set y-axis title
    ),
    xaxis2=dict(
        zeroline=False,
        domain=[0.85, 1],
        showgrid=False,
        showticklabels=False,
    ),
    yaxis2=dict(
        zeroline=False,
        domain=[0.85, 1],
        showgrid=False,
        showticklabels=False,
    ),
    showlegend=False,
    margin=dict(l=0, r=0, t=0, b=0),
    plot_bgcolor='rgba(0,0,0,0)',
    paper_bgcolor='rgba(0,0,0,0)'
)

fig.show()

# Set font for vectorized export
pio.kaleido.scope.font = "Arial"

# Save the figure as a PNG image
pio.write_image(fig, 'test.pdf')



#JUST R VALUE

sns.distplot(Res['R_value'], kde=True, vertical=True, rug=True, hist=False, kde_kws=dict(shade=True), rug_kws=dict(lw=2, color='orange'))
plt.savefig("Density histogram R.pdf")

sns.jointplot(x=Res['Resolution'], y=Res['R_value'], kind='kde',color='#4CB391', xlim=(0,4))
plt.savefig("R Jointplot.pdf")

#JUST FREE R

sns.distplot(Res['R_free'], kde=True, vertical=True, rug=True, hist=False, kde_kws=dict(shade=True), rug_kws=dict(lw=2, color='orange'))
plt.savefig("Density histogram R free.pdf")

sns.jointplot(x=Res['Resolution'], y=Res['R_free'], kind='hex',color='#4CB391', xlim=(0,4))
plt.savefig("R free Jointplot.pdf")



