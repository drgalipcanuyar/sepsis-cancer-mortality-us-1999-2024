from pathlib import Path
import numpy as np
import pandas as pd
from common import *

d = build_primary_strata()
w = pooled_weights(d)

# Save analysis-ready grouped data and standardization weights.
d.to_csv(OUT/'primary_strata_1999_2024.csv', index=False)
w.rename('weight').reset_index().to_csv(OUT/'pooled_age_sex_weights.csv', index=False)

# Direct standardized annual estimates.
c = standardized_series(d, 'cancer', w); c['group']='Cancer-attributed'
n = standardized_series(d, 'noncancer', w); n['group']='Noncancer-attributed'
annual = pd.concat([c,n], ignore_index=True)
annual.to_csv(OUT/'annual_standardized_proportions.csv', index=False)

summary = direct_endpoint_summary(d,w)
rows=[]
for g in ['cancer','noncancer']:
    s=summary[g]
    rows.append([g,s['p1999'],s['p2024'],s['change'],s['change_l'],s['change_u'],s['relative'],s['rel_l'],s['rel_u']])
pd.DataFrame(rows,columns=['group','p1999','p2024','absolute_change','change_lcl','change_ucl','relative_change','relative_lcl','relative_ucl']).to_csv(OUT/'direct_endpoint_summary.csv',index=False)

# Primary grouped-binomial spline model.
m=primary_long(d)
fit,X,cov,disp=fit_spline_model(m,4,False)
wald,df,p=interaction_wald(fit,X,cov)
model_rows=[]
grads={}
for year in range(1999,2025):
    for status in ['cancer','noncancer']:
        pred,se,grad=marginal_prediction(fit,X,cov,w,year,status)
        model_rows.append([year,status,pred,se,pred-1.96*se,pred+1.96*se])
        grads[(year,status)] = grad
model=pd.DataFrame(model_rows,columns=['year','group','proportion','se','lcl','ucl'])
model.to_csv(OUT/'model_marginal_predictions.csv',index=False)
# Differential change 2024 vs 1999.
g = grads[(2024,'cancer')]-grads[(1999,'cancer')]-grads[(2024,'noncancer')]+grads[(1999,'noncancer')]
val=(model.query("year==2024 and group=='cancer'").proportion.iloc[0]-model.query("year==1999 and group=='cancer'").proportion.iloc[0]
     -model.query("year==2024 and group=='noncancer'").proportion.iloc[0]+model.query("year==1999 and group=='noncancer'").proportion.iloc[0])
se=float(np.sqrt(g@cov@g))
with open(OUT/'primary_model_summary.txt','w') as f:
    f.write(f'Grouped observations: {len(m)}\n')
    f.write(f'Model parameters: {X.shape[1]}\n')
    f.write(f'Residual df: {fit.df_resid:.0f}\n')
    f.write(f'Pearson chi-square: {fit.pearson_chi2:.2f}\n')
    f.write(f'Pearson dispersion: {disp:.3f}\n')
    f.write(f'Residual deviance: {fit.deviance:.2f}\n')
    f.write(f'Deviance/df: {fit.deviance/fit.df_resid:.3f}\n')
    f.write(f'Cancer status x time Wald chi-square({df}): {wald:.2f}; P={p:.6g}\n')
    f.write(f'Model differential change 2024 vs 1999: {100*val:.3f} pp (95% CI {100*(val-1.96*se):.3f}, {100*(val+1.96*se):.3f})\n')

print('Primary analysis complete.')
print(f'Interaction Wald chi-square({df}) = {wald:.2f}, P={p:.3g}; dispersion={disp:.3f}')
print(f'Differential model-based change = {100*val:.2f} pp (95% CI {100*(val-1.96*se):.2f} to {100*(val+1.96*se):.2f})')
print('Direct standardized 1999/2024 (%):')
for g in ['cancer','noncancer']:
    s=summary[g]; print(g, round(100*s['p1999'],2), round(100*s['p2024'],2))
