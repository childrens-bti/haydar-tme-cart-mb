# KDM6B remodeling analysis

The canonical KDM6B remodeling workflow is maintained in the
[`haydar-ad-hoc` submodule](../../external/haydar-ad-hoc) at
`analyses/kdm6b-remodeling`.  It replaces the historical KDM6B notebooks in
`analyses/downstream-analyses/` as the source of truth for KDM6B clustering,
abundance, transcriptional, correlation, and immune-program analyses.

Initialize the pinned source repository after cloning this project:

```bash
git submodule update --init --recursive
```

Run the canonical workflow from this repository root inside the project
container:

```bash
bash analyses/kdm6b-remodeling/run-kdm6b-remodeling.sh
```

The launcher changes into the submodule before rendering. This preserves the
source workflow's repository-relative `data/v3`, `figures`, `plots`, and
`results` paths. Generated files remain in the submodule and are not written
to the historical downstream-analysis output directories.
