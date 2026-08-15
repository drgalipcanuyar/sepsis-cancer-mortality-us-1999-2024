import numpy as np
import pandas as pd
from common import *

d=build_primary_strata(); w=pooled_weights(d); m=primary_long(d)
rows=[]
for label,df_t,linear in [('3-df spline',3,False),('Primary 4-df spline',4,False),('5-df spline',5,False),('Linear time',1,True)]:
    fit,X,cov,disp=fit_spline_model(m,df_t,linear)
    stat,df,p=interaction_wald(fit,X,cov)
    vals={}; grads={}
    for y in [1999,2024]:
        for g in ['cancer','noncancer']:
            vals[(y,g)],_,grads[(y,g)]=marginal_prediction(fit,X,cov,w,y,g)
    delta=(vals[(2024,'cancer')]-vals[(1999,'cancer')])-(vals[(2024,'noncancer')]-vals[(1999,'noncancer')])
    grad=grads[(2024,'cancer')]-grads[(1999,'cancer')]-grads[(2024,'noncancer')]+grads[(1999,'noncancer')]
    se=np.sqrt(grad@cov@grad)
    rows.append([label,stat,df,p,disp,100*vals[(1999,'cancer')],100*vals[(2024,'cancer')],100*vals[(1999,'noncancer')],100*vals[(2024,'noncancer')],100*delta,100*(delta-1.96*se),100*(delta+1.96*se)])
pd.DataFrame(rows,columns=['specification','wald_chi2','df','p_value','pearson_dispersion','cancer_1999_pct','cancer_2024_pct','noncancer_1999_pct','noncancer_2024_pct','differential_change_pp','lcl','ucl']).to_csv(OUT/'temporal_specification_sensitivity.csv',index=False)

# Alternative direct-standardization references.
alt=[]
for label,yr in [('Pooled 1999-2024',None),('1999',1999),('2010',2010),('2019',2019)]:
    ww=pooled_weights(d,yr); s=direct_endpoint_summary(d,ww)
    alt.append([label,100*s['cancer']['p1999'],100*s['cancer']['p2024'],100*s['cancer']['change'],100*s['noncancer']['p1999'],100*s['noncancer']['p2024'],100*s['noncancer']['change'],100*s['between']['differential'],s['cancer']['p2024']/s['noncancer']['p2024']])
pd.DataFrame(alt,columns=['reference','cancer_1999_pct','cancer_2024_pct','cancer_change_pp','noncancer_1999_pct','noncancer_2024_pct','noncancer_change_pp','differential_change_pp','ratio_2024']).to_csv(OUT/'alternative_standardization.csv',index=False)
print('Sensitivity analyses complete.')
