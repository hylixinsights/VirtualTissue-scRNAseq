"""Plot saved engineering evidence; does not run or select simulations."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
root=Path(__file__).resolve().parents[1];p=root/'reports/product-v04'
d=json.loads((p/'dynamics-validation.json').read_text());m=d['mean_field'];e=json.loads((p/'ensemble.json').read_text())
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':160})
fig,axes=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
x=np.arange(4)
axes[0].errorbar(x-.07,m['mean'],yerr=np.array(m['standard_error'])*4,fmt='o',color='#146c60',label='64-seed mean ± 4 SE')
axes[0].scatter(x+.07,m['fixed_step_expectation'],marker='x',color='#bf6c67',label='Shared step expectation')
axes[0].set(xticks=x,xticklabels=['Memory','Activated','Cytotoxic','Dysfunctional'],ylabel='Live CD8 count at 12 h',title='Common-kernel verification')
axes[0].tick_params(axis='x',rotation=15);axes[0].legend(fontsize=8)
for i,row in enumerate(e['paired_differences']):
 values=row['paired_tumor_differences'];axes[1].scatter(values,np.arange(5)*.035+i-.07,color='#739e92',s=18)
 axes[1].scatter(row['mean'],i,marker='D',color='#173f39',s=28)
axes[1].axvline(0,color='#abb7b0',lw=1)
axes[1].set(yticks=range(7),yticklabels=['Base','Reactivity ×0.5','Reactivity ×2*','Checkpoint ×0.5','Checkpoint ×2*','Transitions ×0.5','Transitions ×2'],xlabel='Tumor count difference (anti-PD-1 − zero exposure)',title='Assumption sensitivity · five paired seeds')
fig.suptitle('Exploratory implementation verification — no patient validation',fontsize=13)
fig.supxlabel('* Fractions capped at 1. Points are simulation seeds; diamonds are means. Not clinical confidence intervals.',fontsize=9)
fig.savefig(p/'verification.png');fig.savefig(p/'verification.svg');plt.close(fig)
