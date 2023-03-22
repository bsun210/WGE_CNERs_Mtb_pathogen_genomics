#Balaji Sundararaman
#12-22-2022
#usage python3 <multi_sample_genomeCov.tsv> <file_name_prefix>

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys, csv, time

t = pd.read_csv(sys.argv[1], header=None, names=['genome','pos', 10000, 1000, 100,], index_col=0, sep="\t")

#print(t.head())

fig, ax = plt.subplots()

bins = list(np.arange(0,31,1))
bins.append(1000)

#print(bins)

plt.figure(figsize=(5,4))
ax=plt.axes([0.18,0.14,0.78,0.78])

colors=['olivedrab','dodgerblue','crimson',]

ax.hist([t[100],t[1000],t[10000]], bins=bins, histtype="stepfilled", alpha=0.5, label=['100 copy','1000 copy','10000 copy'], density=True, color=colors, cumulative=-1)

xlabels = ['0','5','10','15','20','25','30+']
ax.set_xlim(-1, 31)
ax.set_xticks([xx+0.5 for xx in range(0,31,5)])
ax.set_xticklabels(xlabels, fontsize=10)
ylabels = ['0.0%','10.0%','20.0%','30.0%','40.0%','50.0%','60.0%','70.0%','80.0%','90.0%','100.0%',]#'101.0%']
ax.set_ylim(0, 1.1)
ax.set_yticks([yy/100.0 for yy in range(0,101,10)])
ax.set_yticklabels(ylabels, fontsize=10)

plt.legend(loc='upper right', ncol=1, fontsize=10, facecolor="white", edgecolor="white")
plt.xlabel('Unique Read Depth',fontsize=12)
plt.ylabel('Genome Bases at ≥X Depth',fontsize=12)
plt.title(str(sys.argv[3]),fontsize=14)
plt.grid(visible=True, which='major', axis='both', color='dimgrey', alpha=0.5, linestyle='--', linewidth=0.5)

fn = str(sys.argv[2])+"_genoCov_hist.png"
plt.savefig(fn, dpi=600)

print('Done plotting!')
