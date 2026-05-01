# Murine Brain Tumor CAR-T scRNA-seq Analysis - endpoint analysis module

## Usage

`bash run_module.sh`

## Folder contents

1. `01-endpoint-samples-analysis.Rmd`: focused analysis on endpoint samples to validate trends observed in earlier timepoints.
2. `plots/`: output plots (empty until the module is run).
3. `results/`: intermediate artifacts (ignored by git).
4. `util/`: helper scripts used by the module.

Approximate runtime: 30 minutes on a machine with 32 cores and 128GB of RAM.

## Analysis module directory structure
```
endpoint-analysis/
├── 01-endpoint-samples-analysis.Rmd
├── 01-endpoint-samples-analysis.html
├── README.md
├── run_module.sh
├── plots
├── results
└── util
    ├── run_doubletfinder.R
    └── run_soupx.R
```