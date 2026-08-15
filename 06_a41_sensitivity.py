from pathlib import Path
import numpy as np
import pandas as pd
from common import *
A41=ROOT/'data'/'raw'/'a41'

def rd(name):
    try: d=pd.read_csv(A41/name,dtype=str,encoding='utf-8-sig')
    except UnicodeDecodeError: d=pd.read_csv(A41/name,dtype=str,encoding='latin-1')
    d['Year']=pd.to_numeric(d['Year'],errors='coerce').astype('Int64')
    d['Deaths']=pd.to_numeric(d['Deaths'].replace({'Suppressed':np.nan,'Missing':np.nan}),errors='coerce')
    return d

def pair(a,b):
    da=rd(a); db=rd(b)
    da=da[da['Sex Code'].isin(['F','M']) & da.Year.notna()].copy(); db=db[db['Sex Code'].isin(['F','M']) & db.Year.notna()].copy()
    da['Year']=da.Year.astype(int); db['Year']=db.Year.astype(int)
    return pd.concat([da[da.Year<=2017],db[db.Year>=2018]],ignore_index=True)

d=build_primary_strata()
ca=pair('cancer + A41 only, 1999–2020.csv','cancer + A41 only, 2018–2024.csv')
aa=pair('all deaths + A41 only, 1999–2020.csv','all deaths + A41 only, 2018–2024.csv')
cc=pair('S11_CANCER_A41_0_14_NUM_a.csv','S11_CANCER_A41_0_14_NUM_b.csv')

def aggregate_num(raw, cancer=False):
    x=raw[raw['Ten-Year Age Groups'].isin(AGE_MAP)].copy(); x['age']=x['Ten-Year Age Groups'].map(AGE_MAP)
    older=x[x.age!='0-14'].groupby(['Year','Sex Code','age'],as_index=False)['Deaths'].sum(min_count=1)
    if cancer:
        young=cc[['Year','Sex Code','Deaths']].copy(); young['age']='0-14'
    else:
        yy=x[x.age=='0-14'].groupby(['Year','Sex Code'],as_index=False)['Deaths'].sum(min_count=1); yy['age']='0-14'; young=yy
    return pd.concat([older,young],ignore_index=True).rename(columns={'Sex Code':'sex','Deaths':'num'})
cn=aggregate_num(ca,True); an=aggregate_num(aa,False)
x=d.merge(cn,on=['Year','sex','age']).merge(an,on=['Year','sex','age'],suffixes=('','_all'))
x['cancer_sepsis']=x['num']; x['all_sepsis']=x['num_all']; x['noncancer_sepsis']=x.all_sepsis-x.cancer_sepsis
w=pooled_weights(x)
s=direct_endpoint_summary(x,w)
# Primary model under A41-only definition.
m=primary_long(x); fit,X,cov,disp=fit_spline_model(m,4,False); stat,df,p=interaction_wald(fit,X,cov)
out=pd.DataFrame([
 ['cancer',100*s['cancer']['p1999'],100*s['cancer']['p2024'],100*s['cancer']['change'],100*s['cancer']['change_l'],100*s['cancer']['change_u']],
 ['noncancer',100*s['noncancer']['p1999'],100*s['noncancer']['p2024'],100*s['noncancer']['change'],100*s['noncancer']['change_l'],100*s['noncancer']['change_u']],
],columns=['group','1999_pct','2024_pct','change_pp','lcl','ucl'])
out.to_csv(OUT/'a41_only_direct.csv',index=False)
with open(OUT/'a41_only_model.txt','w') as f:
    f.write(f'Differential direct change: {100*s["between"]["differential"]:.3f} pp (95% CI {100*s["between"]["diff_l"]:.3f}, {100*s["between"]["diff_u"]:.3f})\n')
    f.write(f'Interaction Wald chi-square({df})={stat:.2f}; P={p:.6g}; dispersion={disp:.3f}\n')
print('A41-only sensitivity complete.')
