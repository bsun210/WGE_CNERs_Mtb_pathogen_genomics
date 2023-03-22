#Balaji Sundararaman
#02-12-22

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.ticker as ticker
import numpy as np
from scipy import stats
import matplotlib.image as mpimg
import sys, csv, math
import pandas as pd
from scipy.stats import pearsonr
import seaborn as sns
#sns.set(style="ticks")

#print(sys.argv[1])
input_ls = list(sys.argv[1].split(','))
sample_ls = list(sys.argv[2].split(','))
print(sample_ls)

normCov = {}

for i in range(len(sample_ls)):
    dx = pd.read_csv(input_ls[i],header=None,index_col=1,sep="\t")
    normCov[sample_ls[i]] = dx[3].div(dx[3].mean()).tolist()

#
normCov_df = pd.DataFrame.from_dict(normCov)
#print(normCov_df.head())

mask = np.triu(np.ones_like(normCov_df.corr(method='spearman')))

fig = plt.figure(figsize=(12, 6))
g = fig.add_gridspec(ncols=4, nrows=2,) #gridspec_kw=dict(height_ratios=[4,1], width_ratios=[6,1]))

axd =fig.add_subplot(g[:, 3:])
axa =fig.add_subplot(g[0, 0])
axb =fig.add_subplot(g[0, 1])
axc =fig.add_subplot(g[0, 2])
axe1 =fig.add_subplot(g[1, 0])

axa.set_position([0.05,0.65,0.08,0.3])
axb.set_position([0.18,0.65,0.13,0.3])
axc.set_position([0.36,0.65,0.08,0.3])
axd.set_position([0.55,0.2,0.45,0.85])
axe1.set_position([0.06,0.1,0.33,0.43])

sns.heatmap(normCov_df.corr(method='spearman'), ax=axd, cmap="YlGnBu", fmt ='.2f', annot=True, mask=mask, vmin=0, vmax=1, linewidth=3, linecolor='white', annot_kws={"size": 10}, cbar_kws = dict(use_gridspec=False,location='top',pad=0/1, fraction=0.1,shrink=0.4)) # label='Spearman Rank Correlation',))

sns.scatterplot(data=normCov_df, x=sample_ls[0], y=sample_ls[5], ax=axa,s=2,)
axa.set_xlabel(sample_ls[0], fontsize=11, labelpad=6)
axa.set_ylabel(sample_ls[5], fontsize=12, labelpad=6)
sns.scatterplot(data=normCov_df, x=sample_ls[6], y=sample_ls[11], ax=axb,s=2,)
axb.set_xlabel(sample_ls[6], fontsize=11, labelpad=6)
axb.set_ylabel(sample_ls[11], fontsize=12, labelpad=6)
sns.scatterplot(data=normCov_df, x=sample_ls[0], y=sample_ls[6], ax=axc,s=2,)
axc.set_xlabel(sample_ls[0], fontsize=11, labelpad=6)
axc.set_ylabel(sample_ls[6], fontsize=12, labelpad=6)

cax = axd.figure.axes[-1]
cax.tick_params(labelsize=12)
cax.set_anchor((0.75,-1.5))
cax.set_xlabel('Spearman Rank Correlation', labelpad=10, fontsize=12,)# fontweight='bold')

gc_df=pd.read_csv(sys.argv[3],header=None,index_col=None,sep=",")
print(gc_df.head())
#['GC','Count','H37Rv_WGS','H37Rv_WGE']

axe1.set_xlabel('GC Content', fontsize=12) # fontweight='semibold', fontsize=16, family='Arial', labelpad=15)
#new_y = np.ma.masked_where(gc_df==0, gc_df)
axe1.scatter(gc_df[0],gc_df[2], s=3, color='dodgerblue', label='H37Rv_WGS', zorder=20)
axe1.scatter(gc_df[0],gc_df[3], s=3, color='salmon', label='H37Rv_WGE', zorder=22)
axe1.set_ylabel('Normalized Coverage', fontsize=12) # fontweight='semibold', fontsize=16, family='Arial', labelpad=15)
axe1.axhline(y=1, color='black', alpha = 0.7, lw = 1, linestyle='-.', zorder=1)
axe1.tick_params(axis='both', which='major', labelsize=10, length=6, width=1)
axe1.set_ylim(-0.05, 3.5)
axe1.set_xlim(25,100)
axe1.legend(loc=1, ncol=1, markerscale=2, frameon=False, fontsize=8, markerfirst=False, handletextpad=0.5)

axe2 = axe1.twinx()
axe2.set_position([0.06,0.1,0.33,0.43])
gc_bar = gc_df[1].div(gc_df[1].sum()).tolist()
axe2.bar(gc_df[0], [g*100.0 for g in gc_bar], color='black', alpha=0.35)

gc = np.repeat(gc_df[0], gc_df[1])
#print(gc.quantile([0.25, 0.5, 0.75]), np.mean(gc), np.median(gc))
axe2.axvline(x=gc.quantile(0.05),color='black', alpha = 0.5, lw = 1, linestyle='--', zorder=0)
axe2.text(gc.quantile(0.05)-2.5, 22, '5th %tile GC', fontsize=10, rotation='vertical')
#axe2.axvline(x=gc.quantile(0.5),color='black', alpha = 0.7, lw = 1, linestyle='--')
axe2.axvline(x=gc.quantile(0.95),color='black', alpha = 0.5, lw = 1, linestyle='--', zorder=0)
axe2.text(gc.quantile(0.95)+0.5, 21, '95th %tile GC', fontsize=10, rotation='vertical')
axe2.set_ylim(-0.5, 35)
axe2.set_ylabel('Percent genome \nat given GC', fontsize=12) #fontweight='bold', fontsize=16, family='Sans', labelpad=10)
axe2.tick_params(axis='both', which='major', labelsize=10, length=6, width=1)

for ax in [axa,axb,axc]:
    ax.tick_params(axis = 'both', which = 'major', labelsize=10, length=5, width=1, top=False, right=False, direction='out')
    ax.minorticks_off()
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
#    ax.legend(loc='upper left', ncol=1, fontsize=12, frameon=False)

for ax in [axa,axb,axc,axd,axe1]:
    if ax==axa:
        p=(-0.38,1.14,'A')
        ax.set_xlim(-0.2,3)
        ax.set_ylim(-0.2,6)
        ax.set_xticks([0,1,2,3],['0','1','2','3'])
        ax.set_yticks([0,2,4,6],['0','2','4','6'])
    elif ax==axb:
        p=(-0.25,1.14,'B')
        ax.set_xlim(-0.2,6)
        ax.set_ylim(-0.2,6)
        ax.set_xticks([0,2,4,6],['0','2','4','6'])
        ax.set_yticks([0,2,4,6],['0','2','4','6'])
    elif ax==axc:
        p=(-0.36,1.14,'C')
        ax.set_xlim(-0.1,3)
        ax.set_ylim(-0.1,6)
        ax.set_xticks([0,1,2,3],['0','1','2','3'])
        ax.set_yticks([0,2,4,6],['0','2','4','6'])
    elif ax==axe1:
        p=(-0.13,1.1,'E')
    elif ax==axd:
        p=(-0.18,1.035,'D')
    ax.text(p[0], p[1], p[2], transform=ax.transAxes, fontsize=14, fontweight='bold', va='top', ha='right')

fn = "newFig3_"+str(sys.argv[4])+".png"
plt.savefig(fn, dpi=600)
