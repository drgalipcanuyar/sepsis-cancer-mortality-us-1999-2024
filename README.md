# Sepsis Involvement in Cancer-Attributed Mortality in the United States, 1999–2024

Reproducibility repository for the manuscript:

**Sepsis Involvement in Cancer-Attributed and Noncancer-Attributed Deaths in the United States, 1999–2024**

This repository contains the CDC WONDER aggregate exports and Python code used to reproduce the principal analyses of temporal trends in sepsis involvement among cancer-attributed and noncancer-attributed deaths in the United States.

## Study question

The study evaluates the **proportion of deaths involving sepsis**, rather than population mortality rates or individual-level sepsis incidence.

- **Cancer-attributed deaths:** ICD-10 C00–C97 recorded as the underlying cause of death (UCD).
- **Sepsis involvement:** ICD-10 A40–A41 recorded among multiple causes of death (MCD).
- **Noncancer-attributed deaths:** all deaths minus cancer-attributed deaths within each calendar-year × age-group × sex stratum.
- **Noncancer sepsis-involved deaths:** all sepsis-involved deaths minus cancer-attributed sepsis-involved deaths within the same stratum.

## Data source

Publicly available aggregate mortality data were obtained from the **Centers for Disease Control and Prevention Wide-ranging Online Data for Epidemiologic Research (CDC WONDER) Multiple Cause of Death database**. Data were accessed on **August 13, 2026**.

The final analytic series uses:

- **1999–2017:** final *Multiple Cause of Death, 1999–2020* dataset.
- **2018–2024:** final *Multiple Cause of Death, 2018–2024, Single Race* dataset.
- **2018–2020:** extracted from both sources for overlap validation only and included once in the analytic series.

All CDC WONDER CSV files in the repository root are aggregate exports. No individual-level or identifiable data are included.

## Repository structure

For compatibility with GitHub web upload, this release uses a flat repository layout. Analysis scripts, CDC WONDER aggregate CSV exports, and reproduced output files are stored in the repository root.

Core executable files include:

- `run_all.py`
- `common.py`
- `01_primary_analysis.py` through `07_malignancy_class_analysis.py`
- CDC WONDER aggregate `.csv` exports
- reproduced result `.csv` and `.txt` files

## Primary analytic construction

Primary analyses use nine age categories:

`0–14, 15–24, 25–34, 35–44, 45–54, 55–64, 65–74, 75–84, ≥85 years`

The <1, 1–4, and 5–14 year groups are combined into 0–14 years. Suppressed cells are never assigned zero and are not statistically imputed. Dedicated aggregate pediatric queries are used for suppressed young-age sepsis numerators.

Deaths with age recorded as **Not Stated** are retained in overall annual counts but excluded from age-specific standardization and regression.

Direct age-sex standardization uses the fixed pooled 1999–2024 age-sex distribution of all deaths with stated age as the primary reference distribution.

## Statistical analyses reproduced

The code reproduces the principal manuscript analyses, including:

1. Direct age-sex-standardized annual sepsis-involvement proportions.
2. Primary grouped-binomial logistic regression with a centered 4-df natural cubic spline for calendar time.
3. The global cancer-status × calendar-time interaction with Pearson scale correction for residual overdispersion.
4. Model-based marginal standardized trajectories and differential temporal change.
5. Data-driven segmented log-linear analyses with 0–3 candidate joinpoints selected by BIC.
6. Alternative 3-df, 5-df, and linear calendar-time specifications.
7. Alternative age-sex standardization reference distributions (1999, 2010, and 2019).
8. Age- and sex-specific changes.
9. Pandemic-period sensitivity analyses.
10. Huber–White robust covariance and exchangeable GEE sensitivity analyses.
11. A41.0–A41.9-only sepsis-definition sensitivity analysis.
12. Solid-versus-hematologic malignancy analysis.

Cancer-site-specific raw query outputs and code can be archived as a separate repository release if required; they are not needed to reproduce the primary inference.

## Software

The manuscript analyses used:

- Python 3.13.5
- pandas 2.2.3
- NumPy 2.3.5
- SciPy 1.17.0
- patsy 1.0.2
- statsmodels 0.14.6
- Matplotlib 3.10.8

Install the required packages with:

```bash
python -m pip install -r requirements.txt
```

## Reproducing the analyses

From the repository root, run:

```bash
python run_all.py
```

Reproduced output files are written to the repository root, replacing the included validation copies when rerun.

## Key validation targets

A successful run should reproduce the following primary values (minor differences in the final displayed decimal may occur across numerical environments):

- Direct standardized cancer sepsis involvement: **3.21% in 1999 → 4.69% in 2024**.
- Direct standardized noncancer sepsis involvement: **6.75% → 7.89%**.
- Global cancer-status × calendar-time interaction: **Wald χ²(4) = 525.13; P < .001**.
- Pearson dispersion: **14.402**.
- Model-based cancer marginal estimate: **2.95% → 4.81%**.
- Model-based noncancer marginal estimate: **6.67% → 7.91%**.
- Model-based differential change by 2024: **+0.61 percentage points** (95% CI approximately **0.39 to 0.83**).
- Selected cancer joinpoints: **2013, 2017, 2020**.
- Selected noncancer joinpoints: **2008, 2012, 2015**.
- Solid malignancies: **2.37% → 4.09%**.
- Hematologic malignancies: **9.24% → 10.20%**.
- Malignancy-class × calendar-time interaction: **Wald χ²(4) = 972.69; P < .001**.
- A41-only differential direct change: approximately **+0.29 percentage points**; interaction **Wald χ²(4) = 526.44; P < .001**.

## Important interpretation note

These analyses characterize **sepsis involvement among death certificates attributed to cancer**. They do not estimate sepsis incidence among living patients with cancer, individual-level sepsis risk, or causal effects of cancer treatment, sepsis care, or health-system changes.

## Data availability

The underlying mortality data are publicly available from CDC WONDER. The exact analytic definitions and query construction are described in the manuscript Supplementary Methods and Supplementary Table S1. The raw aggregate query exports used by the scripts are included in this repository to facilitate reproducibility.

## License

Analysis code is released under the MIT License. CDC WONDER data remain subject to the applicable CDC WONDER terms and data-use restrictions.
