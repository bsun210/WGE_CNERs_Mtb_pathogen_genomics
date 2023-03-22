#Balaji Sundararaman
#12-21-22

#Usage
#python3

import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.ticker as ticker
import numpy as np
import sys, csv, math
import pandas as pd
import seaborn as sns
from matplotlib.patches import PathPatch

#dgx = pd.read_csv('MTBC_c10000_drugGenes_normCov.csv', header=0, index_col=3)
#rdx = pd.read_csv('MTBC_c10000_regDiffs_normCov.csv', header=0, index_col=3)
#rdx = pd.read_csv(sys.argv[1], header=0, index_col=3)
rdwgs = pd.read_csv('../B230123_MTBC_WGS/WGS_regDiffs_normCov.csv',header=0, index_col=3, sep=",")
rdwgs.rename({'Mtb_H37Rv':'H37Rv_WGS','Mtb_L1':'L1_WGS','Mtb_L4':'L4_WGS','Mtb_L2':'L2_WGS','Mtb_L3':'L3_WGS','Mbovis':'Mbovis_WGS',}, axis=1, inplace=True)
#print(rdwgs.head())

rdn = pd.read_csv('MTBC_WGE1e4_regDiffs_normCov.csv',header=0, index_col=3, sep=",")
gc = pd.read_csv('MTBC_WGE1e4_picardGC.txt',header=0, index_col=None, sep=",")

#sp = ['M.tb_H37Rv','M.tb_Indo-Oceanic','M.tb_Euro-American','M.tb_East-Asian','M.tb_East-African-Indian','M.bovis',]
sp = ['Mtb_H37Rv','Mtb_L1','Mtb_L4','Mtb_L2','Mtb_L3','Mbovis',]
spw = ['H37Rv_WGS','L1_WGS','L4_WGS','L2_WGS','L3_WGS','Mbovis_WGS']
rdgcn = {s:[] for s in sp}
rdgcn['RD'] = rdwgs.index.tolist()
for s in sp:
    for l in range(77): #len(rdn[s].tolist()):
        rdgcn[s].append(rdn[s][l]/(gc[s][np.around(rdn['gc'][l]*100)]))
#
#print(rdgcn)
#print(bla)

nrd = pd.DataFrame.from_dict(rdgcn)
#nrd.set_index('RD')
nrd.rename({'Mtb_H37Rv':'H37Rv_WGE','Mtb_L1':'L1_WGE','Mtb_L4':'L4_WGE','Mtb_L2':'L2_WGE','Mtb_L3':'L3_WGE','Mbovis':'Mbovis_WGE',}, axis=1, inplace=True)
nrd.set_index('RD',inplace=True)
#print(nrd.head())

rdwgwe = pd.concat([rdwgs.loc[:,spw],nrd], verify_integrity=True, axis=1)#.sort_index().set_index('index')
rdwgwe= rdwgwe.reindex(sorted(rdwgwe.columns), axis=1)
print(rdwgwe.head())

#sns.set_context("talk")
sns.despine()
sns.set_style("ticks")
sns.set_context('talk', rc={'axes.linewidth':1, 'lines.markersize':5, 'axes.labelsize':12, 'xtick.labelsize':10, 'ytick.labelsize':10, 'axes.titlesize':13, 'legend.fontsize':11},)

fig, ax = plt.subplots(figsize=(12, 4))
plt.subplots_adjust(left=0.1, bottom=0.35, right=0.98, top=0.95,)# hspace=0.5)

#g = fig.add_gridspec(ncols=1, nrows=2,)
#axd =fig.add_subplot(g[0, 0])
#axr =fig.add_subplot(g[1, 0])

#sns.heatmap(rdx.loc[:,sp].T,vmin=0,vmax=2, center=1, cmap='bwr', xticklabels=1, yticklabels=1, linewidth=1, cbar_kws={'fraction':0.05, 'pad':0.03})
sns.heatmap(rdwgwe.T,vmin=0,vmax=2, center=1, cmap='bwr', xticklabels=1, yticklabels=1, linewidth=1, cbar_kws={'fraction':0.05, 'pad':0.03})

ax.set_xticklabels(rdwgwe.index.tolist(), fontsize=8)
ax.set_xlabel('Normalized Coverage at Regions of Difference (RD) loci', fontsize=12)

fn = sys.argv[1]+'_RD_HeatMap.png'
plt.savefig(fn, dpi=600)

print('Done plotting!')
