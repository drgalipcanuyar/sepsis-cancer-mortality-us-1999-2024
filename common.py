from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
from scipy.special import expit
from scipy import stats

ROOT = Path(__file__).resolve().parent
CORE = ROOT
OUT = ROOT

AGE9 = ["0-14", "15-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75-84", "85+"]
AGE_MAP = {
    '< 1 year': '0-14', '1-4 years': '0-14', '5-14 years': '0-14',
    '15-24 years': '15-24', '25-34 years': '25-34', '35-44 years': '35-44',
    '45-54 years': '45-54', '55-64 years': '55-64', '65-74 years': '65-74',
    '75-84 years': '75-84', '85+ years': '85+'
}

def read_cdc(path: Path | str) -> pd.DataFrame:
    """Read a CDC WONDER CSV export and coerce Year/Deaths safely."""
    df = pd.read_csv(path, dtype=str, encoding='utf-8-sig')
    df.columns = [c.strip() for c in df.columns]
    if 'Year' in df:
        df['Year'] = pd.to_numeric(df['Year'], errors='coerce').astype('Int64')
    if 'Deaths' in df:
        df['Deaths'] = pd.to_numeric(
            df['Deaths'].replace({'Suppressed': np.nan, 'Missing': np.nan, 'Not Applicable': np.nan}),
            errors='coerce'
        )
    return df

def sex_rows(df: pd.DataFrame) -> pd.DataFrame:
    return df[df['Sex Code'].isin(['F', 'M']) & df['Year'].notna()].copy()

def combine_sources(a: pd.DataFrame, b: pd.DataFrame) -> pd.DataFrame:
    """Use 1999-2017 from legacy final data and 2018-2024 from Single Race final data."""
    a = sex_rows(a)
    b = sex_rows(b)
    a['Year'] = a['Year'].astype(int)
    b['Year'] = b['Year'].astype(int)
    return pd.concat([a[a.Year <= 2017], b[b.Year >= 2018]], ignore_index=True)

def _load_pair(stem_a: str, stem_b: str) -> pd.DataFrame:
    return combine_sources(read_cdc(CORE / stem_a), read_cdc(CORE / stem_b))

def build_primary_strata() -> pd.DataFrame:
    """Construct the 26 x 9 x 2 primary age-sex strata with cancer/noncancer numerators/denominators."""
    A1 = _load_pair('A1a_cancer_sepsis_1999_2020.csv', 'A1b_cancer_sepsis_2018_2024.csv')
    A2 = _load_pair('A2a_all_cancer_deaths_1999_2020.csv', 'A2b_all_cancer_deaths_2018_2024.csv')
    A3 = _load_pair('A3a_all_deaths_with_sepsis_1999_2020.csv', 'A3b_all_deaths_with_sepsis_2018_2024.csv')
    A4 = _load_pair('A4a_all_deaths_1999_2020.csv', 'A4b_all_deaths_2018_2024.csv')
    C1 = _load_pair('C1a_cancer_sepsis_under15_1999_2020.csv', 'C1b_cancer_sepsis_under15_2018_2024.csv')
    C3 = _load_pair('C3a_all_sepsis_under15_1999_2020.csv', 'C3b_all_sepsis_under15_2018_2024.csv')

    def aggregate_den(df: pd.DataFrame, name: str) -> pd.DataFrame:
        d = df[df['Ten-Year Age Groups'].isin(AGE_MAP)].copy()
        d['age'] = d['Ten-Year Age Groups'].map(AGE_MAP)
        out = d.groupby(['Year', 'Sex Code', 'age'], as_index=False)['Deaths'].sum(min_count=1)
        return out.rename(columns={'Deaths': name, 'Sex Code': 'sex'})

    def aggregate_num(df: pd.DataFrame, cdf: pd.DataFrame, name: str) -> pd.DataFrame:
        d = df[df['Ten-Year Age Groups'].isin(AGE_MAP)].copy()
        d['age'] = d['Ten-Year Age Groups'].map(AGE_MAP)
        d = d[d.age != '0-14'].groupby(['Year', 'Sex Code', 'age'], as_index=False)['Deaths'].sum(min_count=1)
        c = cdf[['Year', 'Sex Code', 'Deaths']].copy()
        c['age'] = '0-14'
        out = pd.concat([d, c], ignore_index=True)
        return out.rename(columns={'Deaths': name, 'Sex Code': 'sex'})

    d = aggregate_num(A1, C1, 'cancer_sepsis')
    for x in [aggregate_den(A2, 'cancer_deaths'), aggregate_num(A3, C3, 'all_sepsis'), aggregate_den(A4, 'all_deaths')]:
        d = d.merge(x, on=['Year', 'sex', 'age'], how='inner', validate='one_to_one')
    d['noncancer_sepsis'] = d['all_sepsis'] - d['cancer_sepsis']
    d['noncancer_deaths'] = d['all_deaths'] - d['cancer_deaths']
    d['Year'] = d['Year'].astype(int)
    d['age'] = pd.Categorical(d['age'], categories=AGE9, ordered=True)
    if len(d) != 468 or d.isna().any().any():
        raise ValueError('Primary strata are incomplete; expected 468 complete year-sex-age rows.')
    return d.sort_values(['Year', 'age', 'sex']).reset_index(drop=True)

def pooled_weights(d: pd.DataFrame, reference_year: int | None = None) -> pd.Series:
    q = d if reference_year is None else d[d.Year == reference_year]
    w = q.groupby(['sex', 'age'], observed=True)['all_deaths'].sum()
    return w / w.sum()

def standardized_series(d: pd.DataFrame, group: str, weights: pd.Series | None = None) -> pd.DataFrame:
    if weights is None:
        weights = pooled_weights(d)
    num = 'cancer_sepsis' if group == 'cancer' else 'noncancer_sepsis'
    den = 'cancer_deaths' if group == 'cancer' else 'noncancer_deaths'
    rows = []
    for year in sorted(d.Year.unique()):
        x = d[d.Year == year].set_index(['sex', 'age'])
        p = x[num] / x[den]
        val = float((weights * p).sum())
        var = float((weights.pow(2) * (p * (1-p) / x[den])).sum())
        se = np.sqrt(var)
        rows.append([year, val, se, max(0, val-1.96*se), min(1, val+1.96*se)])
    return pd.DataFrame(rows, columns=['year', 'proportion', 'se', 'lcl', 'ucl'])

def primary_long(d: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in d.iterrows():
        rows.append(dict(year=int(r.Year), age=str(r.age), sex=r.sex, status='cancer', y=r.cancer_sepsis, n=r.cancer_deaths))
        rows.append(dict(year=int(r.Year), age=str(r.age), sex=r.sex, status='noncancer', y=r.noncancer_sepsis, n=r.noncancer_deaths))
    out = pd.DataFrame(rows)
    out['year_c'] = out.year - 1999
    return out

def fit_spline_model(m: pd.DataFrame, df_time: int = 4, linear: bool = False):
    if linear:
        formula = "C(status)*C(age)*C(sex) + year_c + C(status):year_c"
    else:
        spline = f"cr(year_c, df={df_time}, constraints='center')"
        formula = f"C(status)*C(age)*C(sex) + {spline} + C(status):{spline}"
    X = patsy.dmatrix(formula, m, return_type='dataframe')
    endog = np.column_stack([m.y.to_numpy(float), (m.n-m.y).to_numpy(float)])
    fit = sm.GLM(endog, X, family=sm.families.Binomial()).fit()
    dispersion = float(fit.pearson_chi2 / fit.df_resid)
    cov_scaled = fit.cov_params().to_numpy() * dispersion
    return fit, X, cov_scaled, dispersion

def interaction_wald(fit, X: pd.DataFrame, cov_scaled: np.ndarray) -> tuple[float, int, float]:
    if any(':cr(' in c and 'status' in c for c in X.columns):
        idx = [i for i,c in enumerate(X.columns) if ':cr(' in c and 'status' in c]
    else:
        idx = [i for i,c in enumerate(X.columns) if 'C(status)' in c and ':year_c' in c]
    b = fit.params.to_numpy()[idx]
    cov = cov_scaled[np.ix_(idx, idx)]
    stat = float(b @ np.linalg.solve(cov, b))
    df = len(idx)
    p = float(stats.chi2.sf(stat, df))
    return stat, df, p

def marginal_prediction(fit, X: pd.DataFrame, cov_scaled: np.ndarray, weights: pd.Series, year: int, status: str, extra: dict | None = None):
    """Standardized marginal p and delta-method gradient/covariance SE."""
    ps, grads = [], []
    for (sex, age), w in weights.items():
        row = {'year_c': year-1999, 'age': str(age), 'sex': sex, 'status': status}
        if extra:
            row.update(extra)
        new = pd.DataFrame([row])
        xn = patsy.build_design_matrices([X.design_info], new, return_type='dataframe')[0].to_numpy()[0]
        p = float(expit(xn @ fit.params.to_numpy()))
        ps.append(w*p)
        grads.append(w*p*(1-p)*xn)
    pstd = float(np.sum(ps))
    grad = np.sum(np.vstack(grads), axis=0)
    se = float(np.sqrt(grad @ cov_scaled @ grad))
    return pstd, se, grad

def direct_endpoint_summary(d: pd.DataFrame, weights: pd.Series | None = None) -> dict:
    if weights is None:
        weights = pooled_weights(d)
    c = standardized_series(d, 'cancer', weights).set_index('year')
    n = standardized_series(d, 'noncancer', weights).set_index('year')
    out = {}
    for name,s in [('cancer',c),('noncancer',n)]:
        p0,p1 = s.loc[1999,'proportion'], s.loc[2024,'proportion']
        se0,se1 = s.loc[1999,'se'], s.loc[2024,'se']
        change = p1-p0
        se_change = np.sqrt(se0**2+se1**2)
        rr = p1/p0
        se_logrr = np.sqrt((se1/p1)**2+(se0/p0)**2)
        out[name] = dict(p1999=p0,p2024=p1,change=change,change_l=change-1.96*se_change,change_u=change+1.96*se_change,
                         relative=rr-1,rel_l=np.exp(np.log(rr)-1.96*se_logrr)-1,rel_u=np.exp(np.log(rr)+1.96*se_logrr)-1)
    diff = out['cancer']['change']-out['noncancer']['change']
    sed = np.sqrt((c.loc[1999,'se']**2+c.loc[2024,'se']**2)+(n.loc[1999,'se']**2+n.loc[2024,'se']**2))
    rrc = out['cancer']['p2024']/out['cancer']['p1999']; rrn=out['noncancer']['p2024']/out['noncancer']['p1999']
    ratio=rrc/rrn
    selog=np.sqrt((c.loc[2024,'se']/c.loc[2024,'proportion'])**2+(c.loc[1999,'se']/c.loc[1999,'proportion'])**2+
                  (n.loc[2024,'se']/n.loc[2024,'proportion'])**2+(n.loc[1999,'se']/n.loc[1999,'proportion'])**2)
    out['between'] = dict(differential=diff, diff_l=diff-1.96*sed, diff_u=diff+1.96*sed,
                          ratio_relative_changes=ratio, ratio_l=np.exp(np.log(ratio)-1.96*selog), ratio_u=np.exp(np.log(ratio)+1.96*selog))
    return out
