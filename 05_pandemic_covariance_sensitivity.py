import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
from scipy import stats
from common import *

d=build_primary_strata(); w=pooled_weights(d); M=primary_long(d)

def delta_for_model(m,end_year=2024,extra_formula=None):
    if extra_formula is None:
        fit,X,cov,disp=fit_spline_model(m,4,False)
    else:
        formula=extra_formula
        X=patsy.dmatrix(formula,m,return_type='dataframe')
        end=np.column_stack([m.y,m.n-m.y])
        fit=sm.GLM(end,X,family=sm.families.Binomial()).fit()
        disp=fit.pearson_chi2/fit.df_resid; cov=fit.cov_params().to_numpy()*disp
    stat,df,p=interaction_wald(fit,X,cov)
    vals={};gr={}
    for y in [1999,end_year]:
        for g in ['cancer','noncancer']:
            extra = {'pandemic_2020_22': int(2020 <= y <= 2022)} if extra_formula is not None else None
            vals[(y,g)],_,gr[(y,g)]=marginal_prediction(fit,X,cov,w,y,g,extra=extra)
    delta=(vals[(end_year,'cancer')]-vals[(1999,'cancer')])-(vals[(end_year,'noncancer')]-vals[(1999,'noncancer')])
    gg=gr[(end_year,'cancer')]-gr[(1999,'cancer')]-gr[(end_year,'noncancer')]+gr[(1999,'noncancer')]
    se=np.sqrt(gg@cov@gg)
    return stat,df,p,disp,delta,se

rows=[]
for label,mask,endyr in [
    ('Primary 1999-2024',np.ones(len(M),dtype=bool),2024),
    ('Pre-pandemic 1999-2019',M.year<=2019,2019),
    ('Exclude 2020',M.year!=2020,2024),
    ('Exclude 2020-2022',~M.year.isin([2020,2021,2022]),2024),
]:
    stat,df,p,disp,delta,se=delta_for_model(M[mask].copy(),endyr)
    rows.append([label,int(mask.sum()),stat,df,p,disp,100*delta,100*(delta-1.96*se),100*(delta+1.96*se)])

# Pandemic indicator sensitivity: same primary spline plus 2020-2022 indicator and status interaction.
mp=M.copy(); mp['pandemic_2020_22']=mp.year.between(2020,2022).astype(int)
spline="cr(year_c, df=4, constraints='center')"
formula=f"C(status)*C(age)*C(sex) + {spline} + C(status):{spline} + pandemic_2020_22 + C(status):pandemic_2020_22"
stat,df,p,disp,delta,se=delta_for_model(mp,2024,formula)
rows.append(['Full period + 2020-2022 indicator',len(mp),stat,df,p,disp,100*delta,100*(delta-1.96*se),100*(delta+1.96*se)])

pd.DataFrame(rows,columns=['analysis','grouped_n','interaction_chi2','df','p_value','pearson_dispersion','differential_change_pp','lcl','ucl']).to_csv(OUT/'pandemic_sensitivity.csv',index=False)

# Huber-White robust covariance.
spline="cr(year_c, df=4, constraints='center')"; formula=f"C(status)*C(age)*C(sex)+{spline}+C(status):{spline}"
X=patsy.dmatrix(formula,M,return_type='dataframe'); end=np.column_stack([M.y,M.n-M.y])
hc=sm.GLM(end,X,family=sm.families.Binomial()).fit(cov_type='HC0')
idx=[i for i,c in enumerate(X.columns) if ':cr(' in c and 'status' in c]
b=hc.params.to_numpy()[idx]; cv=hc.cov_params().to_numpy()[np.ix_(idx,idx)]; hcstat=float(b@np.linalg.solve(cv,b))

# GEE with 36 status x age x sex clusters, exchangeable working correlation and robust covariance.
mg=M.copy(); mg['cluster']=mg.status+'|'+mg.age+'|'+mg.sex
g=sm.GEE(mg.y/mg.n,X,groups=mg.cluster,family=sm.families.Binomial(),cov_struct=sm.cov_struct.Exchangeable(),weights=mg.n).fit()
b=g.params.to_numpy()[idx]; cv=g.cov_params().to_numpy()[np.ix_(idx,idx)]; geestat=float(b@np.linalg.solve(cv,b))
pd.DataFrame([['Huber-White HC0',hcstat,4,stats.chi2.sf(hcstat,4)],['GEE exchangeable',geestat,4,stats.chi2.sf(geestat,4)]],columns=['variance_specification','interaction_chi2','df','p_value']).to_csv(OUT/'robust_gee_sensitivity.csv',index=False)
print('Pandemic and covariance sensitivity analyses complete.')
