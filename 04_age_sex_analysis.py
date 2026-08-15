import pandas as pd
from common import *
d=build_primary_strata()
rows=[]
for age in AGE9:
    for sex in ['F','M']:
        x=d[(d.age==age)&(d.sex==sex)].set_index('Year')
        for group,num,den in [('cancer','cancer_sepsis','cancer_deaths'),('noncancer','noncancer_sepsis','noncancer_deaths')]:
            p99=x.loc[1999,num]/x.loc[1999,den]; p24=x.loc[2024,num]/x.loc[2024,den]
            rows.append([group,age,sex,100*p99,100*p24,100*(p24-p99)])
pd.DataFrame(rows,columns=['group','age_group','sex','1999_pct','2024_pct','absolute_change_pp']).to_csv(OUT/'age_sex_changes.csv',index=False)
print('Age-sex analysis complete.')
