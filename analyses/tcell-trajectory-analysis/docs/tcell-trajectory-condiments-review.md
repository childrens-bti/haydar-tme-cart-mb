# Condiments outputs

The Rmd renders one compartment at a time (`compartment` = `CD4-like` or `CD8-like`); the module runner renders both. Below, `{cell}` means `cd4` or `cd8`. Results are in `results/cd4-cd8-condiments/`; figures are in `plots/cd4-cd8-condiments/`.

## What is tested, and what is descriptive?

**Condiments tests are cell-level sensitivity analyses.** They test differences among cells from different conditions. They do not provide mouse-level p-values: there are two mice per condition, and mouse identity is not modeled as a replicate effect in these tests. Their p-values and FDRs should not be interpreted as replicate-level evidence.

- **Topology (global, one test per compartment):** Slingshot is refitted on all tested cells (tumor-only cells excluded). Each construct's curves are then refitted on its own cells with the lineage structure held fixed, and compared with curves fitted after 100 permutations of construct labels, using the condiments classifier method (the default for more than two conditions). The statistic is the classifier's accuracy at telling the real-label fits from the permuted ones. A small p-value means at least one construct's curves differ; it does not identify which construct or pair. Because the lineage structure is fixed, the test detects differences in where curves run, not graph-edge changes or branch gain/loss; a nonsignificant result does not establish equivalence.
- **Progression:** tests shared pseudotime distributions globally and on each lineage, using normalized lineage weights.
- **Fate selection:** tests differences in normalized lineage-weight allocations globally and between lineage pairs. It does not establish irreversible terminal fate.

**Mouse-level summaries are descriptive, not additional statistical tests.** Each mouse is shown separately to check direction and consistency. D1 and D2 are mouse identifiers within a condition, not matched mice across treatments. Shared cells contribute to multiple lineages through soft weights; cells are not assigned to one biological program.

## Figures

| File | Pages and purpose | Analysis level |
| --- | --- | --- |
| `tcell_{cell}_condiments_topology_curves.pdf` | The curves compared by the global topology test, drawn in the SLING space where Slingshot fitted them. 1: each lineage, with every construct's curve and the all-condition curve dashed (SLING dimensions 1–2). 2: each construct, with its cells and its own three curves (dimensions 1–2). 3–4: the same two views in dimensions 3–4, where lineages that overlap in dimensions 1–2 can separate. | Cell-level descriptive (no p-values on the figure) |
| `tcell_{cell}_condiments_progression.pdf` | 1: focused condition-pair weighted cumulative distributions on each lineage, with cached KS FDRs and descriptive pooled cell-weighted means. 2: all-comparison progression FDRs. 3: separate bars for each mouse's lineage-weighted mean pseudotime. | Page 1: cell-level descriptive distributions and cell-level sensitivity FDRs; page 2: cell-level sensitivity; page 3: mouse-level descriptive |
| `tcell_{cell}_condiments_fate.pdf` | 1: focused condition-pair cumulative distributions of per-cell relative weights for all three lineage pairs, with cached MMD FDRs and descriptive cell means. 2: all-comparison fate-selection FDRs. 3: each mouse's mean normalized lineage weights as a stacked bar. | Page 1: cell-level descriptive distributions and cell-level sensitivity FDRs; page 2: cell-level sensitivity; page 3: mouse-level descriptive |

All figures exclude tumor-only cells and use the shared condition labels. L1-L3 retain their lineage IDs; biological annotations are in `config/tcell-lineage-annotation-manifest.tsv` and on the descriptive figure panels. An FDR measures statistical evidence, not biological effect size. Progression shifts describe cell positions, not transition speed. A stacked lineage-weight bar is not a count of cells assigned to a terminal fate.

The focused pages default to `plot_comparison = "B7H3_vs_41BB-L"`. This plotting-only Rmd parameter can select another comparison ID in the comparison manifest; it does not change which tests are run or their FDR families/cache contexts. All 10 prespecified comparisons remain in the overviews. The broad mouse summaries also retain all five tested conditions.

**Progression example:** "Within CD8 L3 (exhausted/dysfunctional), are B7H3 cells later than 41BB-L cells?" Inspect the L3 distributions and their descriptive weighted means, alongside the L3 KS FDR. A curve shifted right indicates later pseudotime. The two-sided KS test detects any distribution difference, including crossing/shape differences; significance alone does not establish a one-sided shift or a faster transition. Mean lines/numbers on the focused page pool cells using their lineage weights; they are not a mouse-level test.

**Fate example:** "Do B7H3 cells show stronger exhaustion-associated membership relative to cytotoxic membership than 41BB-L cells?" Inspect the CD8 L1-versus-L3 panel: each cell contributes `w3 / (w1 + w3)` (0 favors L1; 1 favors L3). Cells with zero total weight in that pair are excluded. This matches the within-pair renormalization used by the fate-selection test, not a hard lineage assignment or a fraction of cells in a terminal state. Curves show all eligible cells; the MMD test internally balances condition sizes by subsampling. The pair-specific FDR tests a distribution difference, not the direction or a standalone change in L3's share among all three lineages. A nonsignificant test does not support the proposed enrichment, regardless of a small descriptive mean difference.

Focused distributions show means as dashed lines and numerical labels. Read biological annotations in the facet headers. Mouse summaries remain separate to check whether pooled-cell patterns are consistent between biological samples. Use `plots/cd4-cd8-slingshot/` to inspect the original reference geometry.

**Topology curves:** the HTML report gives the global topology result (method, statistic, p-value and the all-condition refit's lineage count) and prints how each refit lineage was matched to the shared L1–L3 by its end cluster. The PDF shows where the curves behind that result run. Curves that follow the dashed all-condition curve suggest a shared path; curves that separate show where constructs differ; a missing curve means that construct lost the lineage. Curves at sparsely populated lineage ends follow a few outlying cells, so separation there is weak evidence of a different path. Each page is a two-dimensional view of a 20-dimensional space, so differences in other dimensions are not shown. The curves are recomputed for the figure and are not cached.

## Result tables and caches

| File | Contents and purpose | Analysis level |
| --- | --- | --- |
| `tcell_{cell}_condiments_results.tsv` | The main results: 81 rows (1 global topology, 40 progression, 40 fate). Columns: compartment, test, lineage and its program, comparison and its set, the two conditions, method, statistic, p-value, BH-adjusted p-value (`p_adj`, within each comparison), number of mice and cells, `status` and `status_note`. The topology row has `comparison = all_conditions`, its `p_adj` equals its p-value (one test), and `status_note` records the all-condition refit's lineage count. | Cell-level sensitivity |
| `tcell_{cell}_sample_pseudotime_summary.tsv` | One row per mouse and lineage: cell count, weight sum and weighted-mean shared pseudotime. | Mouse-level descriptive |
| `tcell_{cell}_sample_lineage_probabilities.tsv` | One row per mouse and lineage: cell count, normalized weight sum and mean lineage membership. The three means sum to one per mouse. | Mouse-level descriptive |
| `tcell_{cell}_condiments_topology_global_objects.rds` | Raw global topology-test result and cache context, keyed by compartment (`all_conditions`). | Cached cell-level sensitivity results |
| `tcell_{cell}_condiments_progression_objects.rds` | Raw progression-test results and cache context, with global and lineage-specific tests. | Cached cell-level sensitivity results |
| `tcell_{cell}_condiments_fate_objects.rds` | Raw fate-selection results and cache context, with global and lineage-pair tests. | Cached cell-level sensitivity results |

The mouse-summary TSVs retain all six conditions of the shared reference, including tumor-only cells. Tumor-only rows are not used in the prespecified tests or the revised figures. Each RDS is a serialized R list containing `context` and `results`; compatible completed comparisons are reused. Failed calls are not cached. Plotting changes do not invalidate these caches.

`05-tcell-cd4-condiments.html` and `05-tcell-cd8-condiments.html` show the lineage labels, the results table, the global topology table, the topology-curve lineage matching, PDF links and session information. They are reports, not another statistical analysis.
