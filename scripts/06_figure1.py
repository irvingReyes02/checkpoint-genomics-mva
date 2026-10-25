"""Figure 1: BUBR1 variant map + pseudokinase-pocket docking (mean ± s.d. across seeds) + ligand efficiency."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7, 'axes.spines.top': False, 'axes.spines.right': False})
s = pd.read_csv('results/scores.tsv', sep='\t', names=['lig', 'E', 'HA', 'LE'])
r = pd.read_csv('results/reps.tsv', sep='\t', names=['lig', 'seed', 'E'])
g = pd.concat([s[['lig', 'E']], r[['lig', 'E']]]).groupby('lig').E.agg(['mean', 'std', 'count']).reset_index().merge(s[['lig', 'LE']])
g = g.sort_values('mean'); g.to_csv('results/docking_summary.tsv', sep='\t', index=False)
cls = {k: 'flav' for k in ('luteolin', 'apigenin', 'fisetin', 'quercetin', 'resveratrol')} | {'ATP': 'ref', 'staurosporine': 'ref'}
col = {'flav': '#2a9d8f', 'ref': '#7f7f7f', 'tki': '#3d5a80'}
c = [col[cls.get(l, 'tki')] for l in g.lig]
fig = plt.figure(figsize=(7.2, 3.8))
ax = fig.add_axes([0.05, 0.70, 0.9, 0.25]); ax.set_xlim(0, 1050); ax.set_ylim(0, 1); ax.axis('off')
ax.add_patch(plt.Rectangle((0, 0.4), 1050, 0.2, color='#e0e0e0'))
for s0, e0, lab, cc in [(62, 226, 'BUB1 N-term (TPR)', '#bcd4e6'), (766, 1050, 'Pseudokinase domain', '#f4a261')]:
    ax.add_patch(plt.Rectangle((s0, 0.35), e0 - s0, 0.3, color=cc)); ax.text((s0 + e0) / 2, 0.5, lab, ha='center', va='center', fontsize=6.5)
for p, lab in [(737, 'p.Leu737Ter\n(ClinVar P/LP)'), (1002, 'p.Asn1002Lys\n(absent gnomAD; AM 0.92)')]:
    ax.plot([p, p], [0.65, 0.85], color='#c1121f', lw=1); ax.text(p, 0.88, lab, ha='center', va='bottom', fontsize=6, color='#c1121f')
ax.text(0, 0.15, '1', fontsize=6); ax.text(1050, 0.15, '1050 aa', fontsize=6, ha='right'); ax.text(0, 1.1, 'a', fontsize=9, fontweight='bold', transform=ax.transAxes)
ax = fig.add_axes([0.08, 0.20, 0.6, 0.40]); x = np.arange(len(g))
ax.bar(x, g['mean'], yerr=g['std'].fillna(0), color=c, width=0.75, error_kw={'lw': 0.6})
atp = g.loc[g.lig == 'ATP', 'mean'].item(); ax.axhline(atp, ls='--', lw=0.6, color='k'); ax.text(len(g) - 0.5, atp - 0.1, 'ATP', fontsize=6, ha='right', va='top')
ax.set_xticks(x); ax.set_xticklabels(g.lig, rotation=70, ha='right', fontsize=5.5); ax.set_ylabel('Vina score (kcal/mol)'); ax.set_ylim(-8, 0)
ax.text(-0.1, 1.05, 'b', fontsize=9, fontweight='bold', transform=ax.transAxes)
ax = fig.add_axes([0.77, 0.20, 0.21, 0.40]); ax.scatter(g['mean'], -g['LE'], c=c, s=10)
for _, q in g.iterrows():
    if q.lig in ('luteolin', 'fisetin', 'quercetin', 'imatinib', 'ATP'): ax.text(q['mean'], -q['LE'] + 0.008, q.lig, fontsize=5, ha='center')
ax.set_xlabel('Vina score (kcal/mol)'); ax.set_ylabel('Ligand efficiency (-kcal/mol/HA)')
ax.text(-0.35, 1.05, 'c', fontsize=9, fontweight='bold', transform=ax.transAxes)
fig.savefig('figures/Fig1.png', dpi=300); fig.savefig('figures/Fig1.pdf')
