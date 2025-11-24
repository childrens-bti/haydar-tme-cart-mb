# release notes

## current release (v2)
- Data release date: 2025-11-14
- status: available

Added files:
- 10X cellranger raw barcode matrix files for each sample, also reorganized directory structure to separate filtered and raw data.
  - `barcodes.tsv.gz`
  - `features.tsv.gz`
  - `matrix.mtx.gz`
  
```
└── v2
    ├── 41BB-L_pos_T_D1
    │   ├── filtered_bc
    │   │   ├── barcodes.tsv.gz
    │   │   ├── features.tsv.gz
    │   │   └── matrix.mtx.gz
    │   └── raw_bc
    │       ├── barcodes.tsv.gz
    │       ├── features.tsv.gz
    │       └── matrix.mtx.gz
    ├── 41BB-L_pos_T_D2
    │   ├── filtered_bc
    │   │   ├── barcodes.tsv.gz
    │   │   ├── features.tsv.gz
    │   │   └── matrix.mtx.gz
    │   └── raw_bc
    │       ├── barcodes.tsv.gz
    │       ├── features.tsv.gz
    │       └── matrix.mtx.gz
    ├── B7H3_at_endpoint
    │   ├── filtered_bc
    │   │   ├── barcodes.tsv.gz
    │   │   ├── features.tsv.gz
    │   │   └── matrix.mtx.gz
    │   └── raw_bc
    │       ├── barcodes.tsv.gz
    │       ├── features.tsv.gz
    │       └── matrix.mtx.gz
    ├── B7H3_pos_T_D1
    │   ├── filtered_bc
    │   │   ├── barcodes.tsv.gz
    │   │   ├── features.tsv.gz
    │   │   └── matrix.mtx.gz
    │   └── raw_bc
    │       ├── barcodes.tsv.gz
    │       ├── features.tsv.gz
    │       └── matrix.mtx.gz
    ├── B7H3_pos_T_D2
    │   ├── filtered_bc
    │   │   ├── barcodes.tsv.gz
    │   │   ├── features.tsv.gz
    │   │   └── matrix.mtx.gz
    │   └── raw_bc
    │       ├── barcodes.tsv.gz
    │       ├── features.tsv.gz
    │       └── matrix.mtx.gz
    ├── CD28-41BB_at_endpoint
    │   ├── filtered_bc
    │   │   ├── barcodes.tsv.gz
    │   │   ├── features.tsv.gz
    │   │   └── matrix.mtx.gz
    │   └── raw_bc
    │       ├── barcodes.tsv.gz
    │       ├── features.tsv.gz
    │       └── matrix.mtx.gz
    ├── CD28-41BB_pos_T_D1
    │   ├── filtered_bc
    │   │   ├── barcodes.tsv.gz
    │   │   ├── features.tsv.gz
    │   │   └── matrix.mtx.gz
    │   └── raw_bc
    │       ├── barcodes.tsv.gz
    │       ├── features.tsv.gz
    │       └── matrix.mtx.gz
    ├── CD28-41BB_pos_T_D2
    │   ├── filtered_bc
    │   │   ├── barcodes.tsv.gz
    │   │   ├── features.tsv.gz
    │   │   └── matrix.mtx.gz
    │   └── raw_bc
    │       ├── barcodes.tsv.gz
    │       ├── features.tsv.gz
    │       └── matrix.mtx.gz
    ├── CD8-41BB_pos_T_D1
    │   ├── filtered_bc
    │   │   ├── barcodes.tsv.gz
    │   │   ├── features.tsv.gz
    │   │   └── matrix.mtx.gz
    │   └── raw_bc
    │       ├── barcodes.tsv.gz
    │       ├── features.tsv.gz
    │       └── matrix.mtx.gz
    ├── CD8-41BB_pos_T_D2
    │   ├── filtered_bc
    │   │   ├── barcodes.tsv.gz
    │   │   ├── features.tsv.gz
    │   │   └── matrix.mtx.gz
    │   └── raw_bc
    │       ├── barcodes.tsv.gz
    │       ├── features.tsv.gz
    │       └── matrix.mtx.gz
    ├── STOP_pos_T_D1
    │   ├── filtered_bc
    │   │   ├── barcodes.tsv.gz
    │   │   ├── features.tsv.gz
    │   │   └── matrix.mtx.gz
    │   └── raw_bc
    │       ├── barcodes.tsv.gz
    │       ├── features.tsv.gz
    │       └── matrix.mtx.gz
    ├── STOP_pos_T_D2
    │   ├── filtered_bc
    │   │   ├── barcodes.tsv.gz
    │   │   ├── features.tsv.gz
    │   │   └── matrix.mtx.gz
    │   └── raw_bc
    │       ├── barcodes.tsv.gz
    │       ├── features.tsv.gz
    │       └── matrix.mtx.gz
    ├── Tumor_only_D1
    │   ├── filtered_bc
    │   │   ├── barcodes.tsv.gz
    │   │   ├── features.tsv.gz
    │   │   └── matrix.mtx.gz
    │   └── raw_bc
    │       ├── barcodes.tsv.gz
    │       ├── features.tsv.gz
    │       └── matrix.mtx.gz
    └── Tumor_only_D2
    │   ├── filtered_bc
    │   │   ├── barcodes.tsv.gz
    │   │   ├── features.tsv.gz
    │   │   └── matrix.mtx.gz
    │   └── raw_bc
    │       ├── barcodes.tsv.gz
    │       ├── features.tsv.gz
    │       └── matrix.mtx.gz
    └── release-notes.md
```

## archived release (v1)
- Data release date: 2025-11-12
- status: available

Added files:
- `Tumor_only_D1`; 10X v2 data from Non-treated tumor sample 1 (day 20 - BC013)
- `Tumor_only_D2`; 10X v2 data from Non-treated tumor sample 2 (day 20 - BC014)
- `STOP_pos_T_D1`; 10X v2 data from STOP (control) CAR- treated tumor sample 1 (day 20 - BC004)
- `STOP_pos_T_D2`; 10X v2 data from STOP (control) CAR- treated tumor sample 2 (day 20 - BC009)
- `B7H3_pos_T_D1`; 10X v2 data from B7H3 CAR- treated tumor sample 1 (day 20 - BC003)
- `B7H3_pos_T_D2`; 10X v2 data from B7H3 CAR- treated tumor sample 2 (day 20 - BC008)
- `B7H3_at_endpoint`; 10X v2 data from B7H3 CAR- treated tumor sample (at endpoint – BC001)
- `41BB-L_pos_T_D1`; 10X v2 data from 41BBL-B7H3 CAR- treated tumor sample 1 (day 20 - BC005)
- `41BB-L_pos_T_D2`; 10X v2 data from 41BBL-B7H3 CAR- treated tumor sample 2 (day 20 - BC010)
- `CD8-41BB_pos_T_D1`; 10X v2 data from CD8.41BB B7H3 CAR- treated tumor sample 1 (day 20 - BC006)
- `CD8-41BB_pos_T_D2`; 10X v2 data from CD8.41BB B7H3 CAR- treated tumor sample 2 (day 20 - BC011)
- `CD28-41BB_pos_T_D1`; 10X v2 data from CD28.41BB B7H3 CAR- treated tumor sample 1 (day 20 - BC007)
- `CD28-41BB_pos_T_D2`; 10X v2 data from CD28.41BB B7H3 CAR- treated tumor sample 2 (day 20 - BC012)
- `CD28-41BB_at_endpoint`; 10X v2 data from CD28.41BB B7H3 CAR- treated tumor sample (at endpoint – BC002)

```
└── v1
    ├── Tumor_only_D1
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    ├── Tumor_only_D2
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    ├── STOP_pos_T_D1
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    ├── STOP_pos_T_D2
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    ├── B7H3_pos_T_D1
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    ├── B7H3_pos_T_D2
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    ├── B7H3_at_endpoint
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    ├── 41BB-L_pos_T_D1
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    ├── 41BB-L_pos_T_D2
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    ├── CD8-41BB_pos_T_D1
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    ├── CD8-41BB_pos_T_D2
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    ├── CD28-41BB_pos_T_D1
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    ├── CD28-41BB_pos_T_D2
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    ├── CD28-41BB_at_endpoint
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    └── release-notes.md
```