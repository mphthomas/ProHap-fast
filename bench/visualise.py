import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

res = pd.read_csv('bench/results.tsv', sep='\t')
keys = ['date', 'host', 'dataset', 'chrom', 'samples', 'transcripts', 'records', 'arm', 'threads']
tot = (res.groupby(keys, dropna=False, as_index=False)
          .agg(wall=('wall', 'sum'), cpu=('cpu', 'sum'),
               tree_mb=('tree_mb', 'max'), max_mb=('max_mb', 'max')))

fig, ax = plt.subplots(1, 2, figsize=(10, 3))
sns.despine(ax=ax[0])
sns.despine(ax=ax[1])
sns.lineplot(
    data=res[res['step']=='prohap'], y='cpu', x='samples', hue='arm', 
    palette=['#fb3640', '#0a2463'],
    marker='o', markeredgecolor=None, markersize=4,
    ax=ax[0])
ax[0].set_xlabel('n Samples')
ax[0].set_ylabel('CPU Runtime (s)')
ax[0].set_yscale('log')
ax[0].set_xscale('log')
ax[0].legend()
ax[0].grid(True, which='both', linestyle='-', linewidth=0.5, color='lightgray')

sns.lineplot(
    data=res[res['step']=='prohap'], y='max_mb', x='samples', hue='arm', 
    palette=['#fb3640', '#0a2463'],
    marker='o', markeredgecolor=None, markersize=4,
    ax=ax[1])
ax[1].set_xlabel('n Samples')
ax[1].set_ylabel('Memory (MB)')
ax[1].set_yscale('log')
ax[1].set_xscale('log')
ax[1].legend()
ax[1].grid(True, which='both', linestyle='-', linewidth=0.5, color='lightgray')
plt.tight_layout()
plt.savefig('bench/results_cpu_memory.png', dpi=300)
plt.show()
