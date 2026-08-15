# Reproducibility validation

The repository was executed in the preparation environment using Python 3.13.5 and statsmodels 0.14.6.

Verified outputs include:

- Primary grouped observations: 936
- Model parameters: 44
- Residual df: 892
- Pearson chi-square: 12,846.59
- Pearson dispersion: 14.402
- Residual deviance: 12,782.00
- Deviance/df: 14.330
- Cancer-status × calendar-time interaction: Wald chi-square(4) = 525.13
- Model-based differential change, 2024 vs 1999: +0.612 pp (95% CI 0.395 to 0.830)
- Cancer joinpoints: 2013, 2017, 2020
- Noncancer joinpoints: 2008, 2012, 2015
- Huber–White interaction: Wald chi-square(4) = 687.83
- Exchangeable GEE interaction: approximately Wald chi-square(4) = 277.3 (small numerical variation can occur with GEE implementation details)
- A41-only interaction: Wald chi-square(4) = 526.44
- Solid vs hematologic interaction: Wald chi-square(4) = 972.69
