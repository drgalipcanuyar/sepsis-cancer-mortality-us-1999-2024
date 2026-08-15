from pathlib import Path
import subprocess, sys
ROOT=Path(__file__).resolve().parent
scripts=[
    '01_primary_analysis.py',
    '02_segmented_analysis.py',
    '03_sensitivity_analysis.py',
    '04_age_sex_analysis.py',
    '05_pandemic_covariance_sensitivity.py',
    '06_a41_sensitivity.py',
    '07_malignancy_class_analysis.py',
]
for s in scripts:
    print(f'\n=== Running {s} ===',flush=True)
    subprocess.run([sys.executable,str(ROOT/'scripts'/s)],check=True,cwd=ROOT)
print('\nReproducibility analyses completed. See outputs/.')
