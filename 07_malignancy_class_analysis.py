from pathlib import Path
import numpy as np
import pandas as pd
from common import *
BASE=ROOT/'data'/'raw'/'malignancy_class'
AGE8=['0-24','25-34','35-44','45-54','55-64','65-74','75-84','85+']
MAP25={'25-34 years':'25-34','35-44 years':'35-44','45-54 years':'45-54','55-64 years':'55-64','65-74 years':'65-74','75-84 years':'75-84','85+ years':'85+'}

def rd(name):
    for enc in ['utf-8-sig','latin-1']:
        try:
            d=pd.read_csv(BASE/name,dtype=str,encoding=enc); break
        except UnicodeDecodeError: pass
    d['Year']=pd.to_numeric(d['Year'],errors='coerce').astype('Int64')
    d['Deaths']=pd.to_numeric(d['Deaths'].replace({'Suppressed':np.nan,'Missing':np.nan}),errors='coerce')
    return d

def pair(a,b):
    da=rd(a); db=rd(b)
    da=da[da['Sex Code'].isin(['F','M']) & da.Year.notna()].copy(); db=db[db['Sex Code'].isin(['F','M']) & db.Year.notna()].copy()
    da['Year']=da.Year.astype(int); db['Year']=db.Year.astype(int)
    return pd.concat([da[da.Year<=2017],db[db.Year>=2018]],ignore_index=True)

def special(a,b):
    d=pair(a,b)
    return d[['Year','Sex Code','Deaths']].rename(columns={'Sex Code':'sex','Deaths':'value'})

primary=build_primary_strata()
# Total cancer 0-24 denominator for deriving the solid 1999-2017 combined denominator.
tot=primary.copy(); tot['age8']=tot.age.astype(str).replace({'0-14':'0-24','15-24':'0-24'})
tot024=tot[tot.age8=='0-24'].groupby(['Year','sex'],as_index=False).cancer_deaths.sum().rename(columns={'cancer_deaths':'total_cancer_den'})

def build(group):
    tag='SOL' if group=='solid' else 'HEM'
    num=pair(f'S9_{tag}_NUM_a.csv',f'S9_{tag}_NUM_b.csv')
    den=pair(f'S9_{tag}_DEN_a.csv',f'S9_{tag}_DEN_b.csv')
    def older(raw,col):
        x=raw[raw['Ten-Year Age Groups'].isin(MAP25)].copy(); x['age']=x['Ten-Year Age Groups'].map(MAP25)
        return x[['Year','Sex Code','age','Deaths']].rename(columns={'Sex Code':'sex','Deaths':col})
    out=older(num,'y').merge(older(den,'n'),on=['Year','sex','age'])
    # Numerator 0-24 from direct aggregate queries.
    yn=special(f'S9_{tag}_NUM_0_24_a.csv',f'S9_{tag}_NUM_0_24_b.csv').rename(columns={'value':'y'})
    if group=='hematologic':
        nd=special('S9_HEM_DEN_0_24_a.csv','S9_HEM_DEN_0_24_b.csv').rename(columns={'value':'n'})
    else:
        # 2018-2024 direct solid denominator; 1999-2017 = all-cancer denominator - hematologic denominator.
        nb=rd('S9_SOL_DEN_0_24_b.csv'); nb=nb[nb['Sex Code'].isin(['F','M']) & nb.Year.notna()].copy(); nb['Year']=nb.Year.astype(int)
        nb=nb[nb.Year>=2018][['Year','Sex Code','Deaths']].rename(columns={'Sex Code':'sex','Deaths':'n'})
        ha=rd('S9_HEM_DEN_0_24_a.csv'); ha=ha[ha['Sex Code'].isin(['F','M']) & ha.Year.notna()].copy(); ha['Year']=ha.Year.astype(int)
        ha=ha[ha.Year<=2017][['Year','Sex Code','Deaths']].rename(columns={'Sex Code':'sex','Deaths':'hem_n'})
        na=tot024[tot024.Year<=2017].merge(ha,on=['Year','sex']); na['n']=na.total_cancer_den-na.hem_n; na=na[['Year','sex','n']]
        nd=pd.concat([na,nb],ignore_index=True)
    young=yn.merge(nd,on=['Year','sex']); young['age']='0-24'
    out=pd.concat([out,young],ignore_index=True); out['class']=group
    return out

D=pd.concat([build('solid'),build('hematologic')],ignore_index=True)
if D[['y','n']].isna().any().any() or len(D)!=832: raise ValueError(f'Incomplete class strata: n={len(D)}')
# Reference weights from pooled all-death distribution, reaggregated to age8.
r=primary.copy(); r['age8']=r.age.astype(str).replace({'0-14':'0-24','15-24':'0-24'})
w=r.groupby(['sex','age8']).all_deaths.sum(); w=w/w.sum(); w.index.names=['sex','age']
# Direct standardization.
rows=[]
for group in ['solid','hematologic']:
    q=D[D['class']==group]
    for year in [1999,2024]:
        x=q[q.Year==year].set_index(['sex','age']); p=x.y/x.n; val=float((w*p).sum()); var=float((w.pow(2)*p*(1-p)/x.n).sum()); rows.append([group,year,val,np.sqrt(var)])
res=pd.DataFrame(rows,columns=['group','year','proportion','se']); res.to_csv(OUT/'malignancy_class_direct.csv',index=False)
# Model.
m=D.rename(columns={'class':'status'}).copy(); m['year_c']=m.Year-1999; m=m.rename(columns={'Year':'year'})
fit,X,cov,disp=fit_spline_model(m,4,False); stat,df,p=interaction_wald(fit,X,cov)
with open(OUT/'malignancy_class_model.txt','w') as f: f.write(f'Interaction Wald chi-square({df})={stat:.2f}; P={p:.6g}; dispersion={disp:.3f}\n')
print(res.assign(pct=lambda z:100*z.proportion)[['group','year','pct']].to_string(index=False))
print(f'Class x time Wald chi-square({df})={stat:.2f}, P={p:.3g}')
