#Balaji Sundararaman
#02-09-23

#Usage
#python3 mtb_100bp_cov_multiSample.py [per100bp.tsv files] [sample_names] file_name_prefix

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
meanCov = []

for i in range(len(sample_ls)):
    dx = pd.read_csv(input_ls[i],header=None,index_col=1,sep="\t")
    normCov[sample_ls[i]] = dx[3].div(dx[3].mean()).tolist()
    meanCov.append(round(dx[3].mean(),2))

nomrCov_df = pd.DataFrame.from_dict(normCov)

#print(nomrCov_df.head())

def hide_current_axis(*args, **kwds):
    plt.gca().set_visible(False)

def corrfunc(x, y, **kwds):
    cmap = kwds['cmap']
    norm = kwds['norm']
    ax = plt.gca()
    ax.tick_params(bottom=False, top=False, left=False, right=False)
    sns.despine(ax=ax, bottom=True, top=True, left=True, right=True)
    r, _ = pearsonr(x, y)
    facecolor = cmap(norm(r))
    ax.set_facecolor(facecolor)
    lightness = (max(facecolor[:3]) + min(facecolor[:3]) ) / 2
    ax.annotate(f"r={r:.2f}", xy=(.5, .5), xycoords=ax.transAxes,
                color='white' if lightness < 0.7 else 'black', size=26, ha='center', va='center')

g = sns.PairGrid(nomrCov_df, diag_sharey=False)
g.map_lower(plt.scatter, s=10, color='black', alpha=0.5)
g.map_diag(hide_current_axis)#sns.histplot, kde=False, diag_kws={"linewidth": 0, "shade": False})
g.map_upper(corrfunc, cmap=plt.get_cmap('GnBu'), norm=plt.Normalize(vmin=0, vmax=1))
g.fig.subplots_adjust(wspace=0.05, hspace=0.05) 

fn = str(sys.argv[3])+"_mtb_100bp_normCov_corrMatrix.png"
plt.savefig(fn, dpi=600)

