"""Exploratory model ladder with penalized Dirichlet-multinomial likelihood.

This infrastructure is tested on synthetic data. It is not a fitted SPECIAL model.
The caller must supply training-only count arrays; no clinical files are loaded.
"""
from copy import deepcopy
import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.special import gammaln, logsumexp, softmax
from .dynamics import STATES, checked_model, solve

FAMILIES=['M0','M1','M1b','M2','M3']
VARIANTS={'V-T':{'transition'},'V-P':{'division','death'},'V-I':{'influx'},
          'V-TP':{'transition','division','death'},'V-TPI':{'transition','division','death','influx'}}
MEAN_CONTEXT={'contact':1.,'target':1.,'reactive':1.,'tumor_fraction':.5,'myeloid_fraction':.1}


def composition(x):
    x=np.asarray(x,dtype=float)
    if not np.isfinite(x).all() or (x<0).any() or (x.sum(axis=-1)<=0).any():raise ValueError('Positive total abundance required.')
    return x/x.sum(axis=-1,keepdims=True)


def clr(x):
    x=np.log(np.clip(composition(x),1e-12,1));return x-x.mean(axis=-1,keepdims=True)


def dm_logpmf(counts,probabilities,concentration):
    n=np.asarray(counts,dtype=float);p=np.asarray(probabilities,dtype=float)
    if concentration<=0 or not np.isfinite(concentration) or n.shape[-1]!=4 or not np.isfinite(n).all() or (n<0).any() or not np.all(n==np.floor(n)):
        raise ValueError('DM requires four integer counts and positive concentration.')
    if not np.isfinite(p).all() or (p<0).any() or not np.allclose(p.sum(axis=-1),1):raise ValueError('Invalid composition.')
    alpha=np.clip(p,1e-12,1);alpha=alpha/alpha.sum(axis=-1,keepdims=True)*concentration;total=n.sum(axis=-1)
    return gammaln(total+1)-gammaln(n+1).sum(axis=-1)+gammaln(concentration)-gammaln(total+concentration)+(gammaln(n+alpha)-gammaln(alpha)).sum(axis=-1)


def parameter_value(m,name,value=None):
    parts=name.split(':');kind=parts[0]
    if kind=='transition' and len(parts)==3:
        matches=[e for e in m['transitions'] if e['source']==parts[1] and e['target']==parts[2]]
        if len(matches)!=1:raise ValueError('Parameter names must identify one permitted edge.')
        if value is not None:matches[0]['rate_per_hour']=float(value)
        return matches[0]['rate_per_hour']
    mapping={'division':'division_per_hour','death':'death_per_hour','influx':'influx_cells_per_hour'}
    if kind not in mapping or len(parts)!=2 or parts[1] not in STATES:raise ValueError('Unknown parameter name.')
    if value is not None:m[mapping[kind]][parts[1]]=float(value)
    return m[mapping[kind]][parts[1]]


def predict(family, baseline, times, *, model=None, contexts=None, shift=None, mean_time=1.):
    if family not in FAMILIES:raise ValueError('Unknown ladder family.')
    baseline=np.asarray(baseline,dtype=float);times=np.asarray(times,dtype=float)
    if baseline.ndim!=2 or baseline.shape[1]!=4 or times.shape!=(len(baseline),) or (times<0).any():raise ValueError('Patient x four-state baseline and one nonnegative time per patient required.')
    p=composition(baseline)
    if family=='M0':return p
    if family in ['M1','M1b']:
        shift=np.zeros(4) if shift is None else np.asarray(shift)
        scale=np.ones(len(times)) if family=='M1' else times/mean_time
        return softmax(clr(p)+scale[:,None]*shift,axis=-1)
    m=checked_model(model);contexts=contexts or [MEAN_CONTEXT]*len(times)
    return np.array([composition(solve(n,[float(t)],m,context if family=='M3' else MEAN_CONTEXT)[0]) for n,t,context in zip(baseline,times,contexts)])


def fit_ladder(pre,on,times,*,family='M3',model=None,parameters=(),contexts=None,concentration=100.,
               baseline_draws=16,regularization=.01,seed=21,patient_ids=None,forbidden_ids=(),variant='V-TPI',evidence='synthetic'):
    if evidence not in ['synthetic','exploratory_training']:raise ValueError('No clinical-validation claim is supported.')
    pre=np.asarray(pre,dtype=float);on=np.asarray(on,dtype=float);times=np.asarray(times,dtype=float)
    if pre.shape!=on.shape or pre.ndim!=2 or pre.shape[1]!=4 or len(pre)==0 or (pre<0).any() or not np.isfinite(pre).all() or (pre!=np.floor(pre)).any():raise ValueError('Paired patient x 4 count matrices required.')
    dm_logpmf(on,np.full_like(on,.25),concentration)
    if times.shape!=(len(pre),) or not np.isfinite(times).all() or (times<=0).any():raise ValueError('Individual positive time intervals required; do not invent missing intervals.')
    if not isinstance(baseline_draws,int) or not 2<=baseline_draws<=1024 or regularization<0:raise ValueError('Use 2–1024 baseline draws and nonnegative regularization.')
    if patient_ids is None:patient_ids=[f'synthetic_{i}' for i in range(len(pre))] if evidence=='synthetic' else None
    if patient_ids is None or len(patient_ids)!=len(pre) or len(set(patient_ids))!=len(patient_ids):raise ValueError('One unique training patient identity per row is required.')
    if set(patient_ids)&set(forbidden_ids):raise ValueError('Leakage guard: forbidden/test patient in training data.')
    if variant not in VARIANTS:raise ValueError('Unknown mechanism variant.')
    if any(n.split(':')[0] not in VARIANTS[variant] for n in parameters):raise ValueError('Parameter is not permitted by mechanism variant.')
    if len(set(parameters))!=len(parameters):raise ValueError('Duplicate parameters.')
    m=checked_model(model)
    # Explicitly remove disallowed mechanisms in each synthetic variant.
    if 'transition' not in VARIANTS[variant]:m['transitions_enabled']=False
    for kind,key in [('division','division_per_hour'),('death','death_per_hour'),('influx','influx_cells_per_hour')]:
        if kind not in VARIANTS[variant]:m[key]=dict.fromkeys(STATES,0.)
    contexts=contexts or [MEAN_CONTEXT]*len(pre)
    if len(contexts)!=len(pre):raise ValueError('One context per patient required.')
    rng=np.random.default_rng(seed)
    baseline=np.array([rng.dirichlet(n+.5,baseline_draws)*max(1,n.sum()) for n in pre])
    mean_time=float(times.mean())
    if family in ['M1','M1b']:
        empirical=(clr(on+.5)-clr(pre+.5)).mean(axis=0);x0=empirical[:3]-empirical[3];bounds=[(-10,10)]*3
    elif family in ['M2','M3']:
        x0=np.log([max(1e-6,parameter_value(m,n)) for n in parameters]);bounds=[(-14,np.log(1.98 if not n.startswith('influx') else 99.)) for n in parameters]
    elif family=='M0':x0=np.array([]);bounds=[]
    else:raise ValueError('Unknown ladder family.')
    def unpack(theta):
        current=deepcopy(m);shift=np.r_[theta,0.] if family in ['M1','M1b'] else np.zeros(4)
        if family in ['M2','M3']:
            for name,value in zip(parameters,np.exp(theta)):parameter_value(current,name,value)
        return current,shift
    def predictions(theta,draws=True):
        current,shift=unpack(theta)
        if not draws:return predict(family,pre+.5,times,model=current,contexts=contexts,shift=shift,mean_time=mean_time)
        # Linear compartment solution is a linear map: evaluate all baseline draws together.
        out=[]
        for i,t in enumerate(times):
            if family in ['M0','M1','M1b']:
                out.append(predict(family,baseline[i],np.full(baseline_draws,t),shift=shift,mean_time=mean_time))
            else:
                ctx=contexts[i] if family=='M3' else MEAN_CONTEXT
                zero=solve(np.zeros(4),[t],current,ctx)[0]
                matrix=np.column_stack([solve(np.eye(4)[j],[t],current,ctx)[0]-zero for j in range(4)])
                out.append(composition(baseline[i]@matrix.T+zero))
        return np.array(out)
    def objective(theta):
        pp=predictions(theta);ll=dm_logpmf(on[:,None,:],pp,concentration)
        return float(-np.sum(logsumexp(ll,axis=1)-np.log(baseline_draws))+.5*regularization*np.sum((theta-x0)**2))
    fit=minimize(objective,x0,method='L-BFGS-B',bounds=bounds) if len(x0) else None
    theta=fit.x if fit else x0;current,shift=unpack(theta);central=predictions(theta,False)
    jac=np.empty((pre.size,len(theta)));eps=1e-4
    for j in range(len(theta)):
        plus=theta.copy();minus=theta.copy();plus[j]+=eps;minus[j]-=eps
        jac[:,j]=((predictions(plus,False)-predictions(minus,False))/(2*eps)).ravel()
    singular=np.linalg.svd(jac,compute_uv=False) if len(theta) else np.array([])
    rank=int(np.linalg.matrix_rank(jac,tol=1e-6)) if len(theta) else 0
    # Penalized observed Hessian, local approximation only; not an identification proof.
    hessian=np.empty((len(theta),len(theta)));h=1e-3
    for i in range(len(theta)):
        for j in range(len(theta)):
            ei=np.eye(len(theta))[i]*h;ej=np.eye(len(theta))[j]*h
            hessian[i,j]=(objective(theta+ei+ej)-objective(theta+ei-ej)-objective(theta-ei+ej)+objective(theta-ei-ej))/(4*h*h)
    covariance=np.linalg.pinv(hessian) if len(theta) else hessian
    se=np.sqrt(np.maximum(0,np.diag(covariance)))
    estimates={name:{'value':float(np.exp(value)),'evidence_class':'estimated','source':evidence,
                    'unit':'cells/hour' if name.startswith('influx') else '1/hour',
                    'local_95_interval':np.exp(np.clip([value-1.96*sd,value+1.96*sd],-30,30)).tolist()}
               for name,value,sd in zip(parameters,theta,se)} if family in ['M2','M3'] else {}
    current['evidence_class']='estimated';current['source']=f'Penalized likelihood on {evidence} training counts; not validated'
    for name in parameters:
        if name.startswith('transition'):
            for edge in current['transitions']:
                if name==f'transition:{edge["source"]}:{edge["target"]}':edge['evidence_class']='estimated';edge['provenance']=current['source']
    current['parameter_provenance']['fitted_subset']=estimates
    return {'family':family,'variant':variant,'backend':f'SciPy {scipy.__version__} L-BFGS-B penalized maximum likelihood',
      'success':bool(fit.success) if fit else True,'message':str(fit.message) if fit else 'No fitted parameters',
      'objective':objective(theta),'estimates':estimates,'shift':shift.tolist(),'mean_training_time':mean_time,
      'concentration':concentration,'baseline_uncertainty':f'{baseline_draws} fixed Dirichlet(n_pre+0.5) draws per patient',
      'regularization':regularization,'predicted_compositions':central.tolist(),'model':current,
      'identifiability':{'jacobian_rank':rank,'parameters':len(theta),'singular_values':singular.tolist(),
       'full_local_rank':rank==len(theta),'note':'Composition cannot distinguish a common growth shift; birth/death with equal net growth are confounded. Local rank is necessary, not sufficient. Intervals are local Hessian approximations.'},
      'patient_ids':patient_ids,'seed':seed,'evidence':evidence,'clinical_validation':False}
