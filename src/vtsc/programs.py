"""Fixed, exploratory human CD8 rank programs; no patient-fitted reference."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata, spearmanr
from .dynamics import STATES

LIBRARY = {
 'version':'cd8-rank-v1', 'evidence_class':'inferred, provisional marker recipe',
 'source':'A05; Li Cell 2019 doi:10.1016/j.cell.2018.11.043; Sade-Feldman Cell 2018 doi:10.1016/j.cell.2018.10.038',
 'method':'Capped positive-gene within-cell rank enrichment, U-statistic normalization; not calibrated state probabilities',
 'state': {'Memory':['TCF7','CCR7','SELL','LEF1'], 'Activated':['GZMK','CD69','IL2RA','TNFRSF9'],
  'Cytotoxic':['GZMB','PRF1','CTSW','FGFBP2'], 'Dysfunctional':['PDCD1','HAVCR2','LAG3','TOX','ENTPD1']},
 'held_out': {'Memory':['IL7R','BCL2'], 'Activated':['NR4A2','EGR2'],
  'Cytotoxic':['GZMH','CCL5'], 'Dysfunctional':['CXCL13','LAYN','TIGIT']},
 'axes':{'dysfunction':['PDCD1','HAVCR2','LAG3','TOX','ENTPD1'],
  'proliferation':['MKI67','TOP2A','PCNA','MCM2','TYMS'],
  'cytotoxicity':['GZMB','PRF1','CTSW','FGFBP2'], 'persistence':['IL7R','BCL2','TCF7']},
 'controls':{'dissociation':['FOS','JUN','JUNB','HSPA1A','EGR1'], 'tcr_signal':['NR4A1','NR4A2','NR4A3','EGR2'],
  'myeloid_inflammation':['CXCL9','CXCL10','STAT1','IRF1'], 'antigen_presentation':['HLA-A','HLA-B','B2M','TAP1']},
 'max_rank':1500,'min_detected_state_markers':2,'min_margin':.08,'min_state_score':.15,
 'limits':'Not Tpex diagnosis, CNV validation, exhaustion protein measurement, or a training-fold reference.'}


def gene_names(data):
    names=data.var.get('gene_symbols',pd.Series(data.var_names,index=data.var_names)).astype(str).to_numpy()
    if np.mean(np.char.startswith(names.astype(str),'ENSG'))>.5:
        raise ValueError('Ensembl gene IDs require a gene-symbol mapping (--gene-map or var.gene_symbols).')
    return names


def score_programs(data):
    names=gene_names(data); x=sparse.csr_matrix(data.X); n=data.n_obs
    panels={f'{group}:{key}':genes for group in ['state','held_out','axes','controls'] for key,genes in LIBRARY[group].items()}
    coverage={k:[g for g in genes if g in names] for k,genes in panels.items()}
    # For duplicate symbols, use the best positive feature rank; never change raw feature IDs.
    needed=sorted(set(sum(panels.values(),[]))); symbol_index={g:np.flatnonzero(names==g) for g in needed}
    scores={k:np.full(n,np.nan) for k in panels}; detected={k:np.zeros(n,dtype=int) for k in panels}
    rmax=min(LIBRARY['max_rank'],data.n_vars)
    for i in range(n):
        row=x.getrow(i); counts=dict(zip(row.indices,row.data)); ranks=dict(zip(row.indices,rankdata(-row.data,method='average')))
        for k,genes in coverage.items():
            if len(genes)<2: continue
            rs=[]
            for g in genes:
                ids=symbol_index[g]; active=[j for j in ids if counts.get(j,0)>0]
                detected[k][i]+=bool(active)
                rs.append(min(rmax,min((ranks[j] for j in active),default=rmax)))
            size=len(rs); offset=size*(size+1)/2
            scores[k][i]=np.clip(1-(sum(rs)-offset)/max(1,size*rmax-offset),0,1)
    return scores,detected,coverage


def classify(scores,detected,group='state'):
    values=np.column_stack([scores[f'{group}:{s}'] for s in STATES]); safe=np.nan_to_num(values,nan=-1)
    order=np.argsort(safe,axis=1,kind='stable'); rows=np.arange(len(values)); top=order[:,-1]
    margin=safe[rows,top]-safe[rows,order[:,-2]]; support=np.column_stack([detected[f'{group}:{s}'] for s in STATES])[rows,top]
    labels=np.asarray(STATES,dtype=object)[top]
    labels[(support<2)|(margin<LIBRARY['min_margin'])|(safe[rows,top]<LIBRARY['min_state_score'])]='Unknown'
    return labels,margin


def audit_programs(data, lineages, scores, detected, coverage, assigned):
    cd8=np.asarray(lineages)=='CD8_T'; held,margin=classify(scores,detected,'held_out')
    valid=cd8&(np.asarray(assigned)!='Unknown')&(held!='Unknown')
    depth=np.asarray((sparse.csr_matrix(data.X)>0).sum(axis=1)).ravel()
    panels={f'{g}:{k}':set(v) for g in ['state','axes','controls'] for k,v in LIBRARY[g].items()}
    overlap={a:{b:len(ga&gb)/len(ga|gb) for b,gb in panels.items()} for a,ga in panels.items()}
    correlations={}
    for key,values in scores.items():
        ok=cd8&np.isfinite(values)
        correlations[key]={}
        for cov,arr in [('detected_genes',depth),('dissociation',scores['controls:dissociation'])]:
            use=ok&np.isfinite(arr)
            correlations[key][cov]=(float(spearmanr(values[use],arr[use]).statistic)
              if use.sum()>3 and np.ptp(values[use])>0 and np.ptp(arr[use])>0 else None)
    kappa=None
    if valid.any():
        agreement=float(np.mean(np.asarray(assigned)[valid]==held[valid]))
        chance=sum(np.mean(np.asarray(assigned)[valid]==s)*np.mean(held[valid]==s) for s in STATES)
        kappa=(agreement-chance)/(1-chance) if chance<1 else None
    return {'library':LIBRARY,'coverage':coverage,'overlap_jaccard':overlap,'cd8_cells':int(cd8.sum()),
      'cd8_state_counts':{s:int(np.sum(cd8&(np.asarray(assigned)==s))) for s in STATES+['Unknown']},
      'heldout_markers':{'assessable_cells':int(valid.sum()),'concordance':float(np.mean(np.asarray(assigned)[valid]==held[valid])) if valid.any() else None,'kappa':kappa,
       'scope':'Independent markers, not independent biological labels; no thresholds tuned to concordance'},
      'confounder_spearman':correlations,'depth_quantiles':np.quantile(depth,[0,.1,.5,.9,1]).tolist(),
      'library_sha256':hashlib.sha256(json.dumps(LIBRARY,sort_keys=True).encode()).hexdigest()}
