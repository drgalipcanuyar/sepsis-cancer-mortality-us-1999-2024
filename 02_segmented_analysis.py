import itertools
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from common import *

def candidate_configs(nj):
    if nj==0: return [()]
    # At least two observed annual points strictly between a joinpoint and endpoint/adjacent joinpoint.
    return itertools.combinations(range(2002, 2022), nj)

def design(years,jps):
    x=np.asarray(years,float)
    return np.column_stack([np.ones(len(x)),x]+[np.maximum(0,x-j) for j in jps])

def fit_config(df,jps):
    y=np.log(df.proportion.to_numpy())
    varlog=(df.se.to_numpy()/df.proportion.to_numpy())**2
    X=design(df.year,jps)
    fit=sm.WLS(y,X,weights=1/varlog).fit()
    return fit

def valid(jps):
    pts=(1999,)+tuple(jps)+(2024,)
    return all(pts[i+1]-pts[i]>=3 for i in range(len(pts)-1))

def select(df):
    table=[]; best_by_n={}
    for nj in range(4):
        best=None
        for jps in candidate_configs(nj):
            if not valid(jps): continue
            fit=fit_config(df,jps)
            if best is None or fit.bic<best[0]: best=(fit.bic,tuple(jps),fit)
        best_by_n[nj]=best
        table.append([nj,best[0],best[1]])
    selected=min(best_by_n.values(), key=lambda z:z[0])
    return selected,table

def segment_results(jps,fit):
    bounds=(1999,)+jps+(2024,)
    out=[]
    for s in range(len(jps)+1):
        c=np.zeros(len(fit.params)); c[1]=1
        for j in range(s): c[2+j]=1
        b=float(c@fit.params); se=float(np.sqrt(c@fit.cov_params()@c)); t=stats.t.ppf(.975,fit.df_resid)
        lo,hi=b-t*se,b+t*se
        p=2*stats.t.sf(abs(b/se),fit.df_resid)
        out.append([bounds[s],bounds[s+1],100*(np.exp(b)-1),100*(np.exp(lo)-1),100*(np.exp(hi)-1),p])
    return out

d=build_primary_strata(); w=pooled_weights(d)
all_bic=[]; all_seg=[]
for group in ['cancer','noncancer']:
    s=standardized_series(d,group,w)
    sel,table=select(s)
    bic,jps,fit=sel
    for nj,b,j in table: all_bic.append([group,nj,b,','.join(map(str,j))])
    for r in segment_results(jps,fit): all_seg.append([group,*r])
    print(group,'joinpoints',jps,'BIC',round(bic,3))
pd.DataFrame(all_bic,columns=['group','joinpoints_n','BIC','best_joinpoint_years']).to_csv(OUT/'segmented_BIC.csv',index=False)
pd.DataFrame(all_seg,columns=['group','start','end','APC_percent','lcl','ucl','p_value']).to_csv(OUT/'segmented_APC.csv',index=False)
