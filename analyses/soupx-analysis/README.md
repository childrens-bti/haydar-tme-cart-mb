# Determination of best SoupX parameters to run on murine CAR-T brain tumor scRNA-seq samples

## Usage

`bash run_module.sh`

## Folder contents

1. `01-ambient-assessment.Rmd`: Examining the degree of ambient RNA contamination with built in SoupX functions using one sample (Tumor_D1). Determine key lineage genes that we can leverage to aid in decontamination.
2. `02-one-sample-test.Rmd`: Using Tumor_D1 sample as example, test the effects of SoupX with 4 parameter sets:
    - No SoupX
    - SoupX with default settings, with automatic estimation of contamination fraction
    - SoupX with contamination fraction manually set to 20%
    - SoupX with manual estimation of contamination fraction using lineage gene sets
3. `03-all-sample-test.Rmd`: Check before vs after running SoupX with lineage gene sets guidance using all samples.

## Analysis module directory structure
```
soupx_analysis/
├── results/
│   └── soupx_clusters_B7H3_pos_T_D2.csv
│
├── util/
│   ├── run_doubletfinder.R
│   └── run_soupx.R
│
├── 01-ambient-assessment.Rmd
├── 01-ambient-assessment.html
├── 02-one-sample-test.Rmd
├── 02-one-sample-test.html
├── 03-all-sample-test.Rmd
├── 03-all-sample-test.html
├── README.md
└── run_module.sh
```