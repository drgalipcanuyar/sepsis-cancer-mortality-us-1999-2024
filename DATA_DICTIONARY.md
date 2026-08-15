# Data dictionary and query labels

The A/B/C labels are study-specific labels and are not CDC WONDER designations.

| Label | UCD restriction | MCD restriction | Grouping | Role |
|---|---|---|---|---|
| A1 | C00–C97 | A40–A41 | Year × 10-y age × sex | Cancer-attributed deaths involving sepsis |
| A2 | C00–C97 | None | Year × 10-y age × sex | All cancer-attributed deaths |
| A3 | None | A40–A41 | Year × 10-y age × sex | All deaths involving sepsis |
| A4 | None | None | Year × 10-y age × sex | All deaths |
| B1 | C00–C97 | A40–A41 | Year | Exact annual cancer-attributed sepsis-involved totals |
| B2 | C00–C97 | None | Year | Exact annual cancer-attributed death totals |
| B3 | None | A40–A41 | Year | Exact annual all-sepsis totals |
| C1 | C00–C97 | A40–A41 | Year × sex, ages <1/1–4/5–14 selected | Directly queried 0–14 cancer sepsis numerator |
| C3 | None | A40–A41 | Year × sex, ages <1/1–4/5–14 selected | Directly queried 0–14 all-sepsis numerator |

The final series takes 1999–2017 from the 1999–2020 final dataset and 2018–2024 from the 2018–2024 Single Race final dataset. The 2018–2020 overlap is used for validation only.
