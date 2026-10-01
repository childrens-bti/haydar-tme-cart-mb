# CD8-like shared trajectory review

Review date: 2026-09-29

Consistency check: 2026-10-01

Scope: descriptive review of the existing shared `all_conditions` CD8-like results from `03-tcell-cd4-cd8-slingshot.Rmd`. No analysis was rerun.

## Purpose

The trajectory orders CD8-like cells along connected transcriptional states to identify possible memory-like, cytotoxic, proliferating, and dysfunctional branches. Pseudotime is an inferred state ordering, not direct observation of cells changing over time.

## Overview

- 2,138 cells, seven retained clusters, and three lineages.
- Cluster 6 was selected as the root using ProjecTILs naïve/memory evidence.
- Cluster 6 is 79.7% ProjecTILs `CD8_NaiveLike` and 79.0% previously annotated naïve/central-memory.
- Early pseudotime has high `Tcf7`, `Lef1`, `Ccr7`, and `Il7r` expression and high naïve/memory AUCell activity.

The CD8 root is biologically well supported. Across the three lineages, naïve/memory identity generally declines while cytotoxic, proliferative, or dysfunction-associated programs increase.

## Lineages

| Lineage | Main result | Supporting evidence | Main caution |
|---|---|---|---|
| 1: Tex-enriched cytotoxic trajectory | Cluster 6 to cluster 1 | `Gzmb` and `Prf1` rise; `Pdcd1`, `Tox`, `Lag3`, and `Havcr2` peak; endpoint is 59.0% ProjecTILs `CD8_Tex`; median pseudotime increases from 9.4 to 33.2 | Endpoint also has NKT/NK-like SingleR support and may not be a pure exhausted-CD8 population |
| 2: proliferating effector/dysfunction trajectory | Cluster 6 to cluster 5 | `Mki67`, `Top2a`, proliferation AUCell, cytotoxicity, and exhaustion-associated programs increase; median pseudotime increases from 8.8 to 46.6 | Cell cycle may drive much of this branch's separation |
| 3: terminal cytotoxic trajectory with exhaustion-associated activity | Cluster 6 to cluster 3 | `Gzmb` and `Ccl5` peak earlier; `Prf1`, `Ifng`, `Pdcd1`, `Tox`, and `Lag3` rise toward the endpoint; exhaustion AUCell rises late; median pseudotime increases from 8.8 to 27.6 | The curve shares much of its early route with lineage 1, so the two terminal programs are not cleanly separated |

## Condition representation

All six conditions have cells with non-missing pseudotime on all three lineages, but their representation is unequal. B7H3 contributes many cells with lineage 1 pseudotime, while CD8-41BB contributes many cells with lineage 3 pseudotime. Cells can contribute to multiple lineages through Slingshot weights, so these are not exclusive fate counts. The condition-specific smooth curves are descriptive, not evidence of condition-specific fate preference.

## Main limitations

- The proliferating branch may be partly defined by cell cycle.
- Lineage 1 may contain an NK/NKT-like cytotoxic population.
- Lineages 1 and 3 share much of their early route.
- The gene curves are descriptive smoothers without uncertainty intervals or significance tests.
- Unequal cell recovery and only two samples per condition limit condition-level interpretation.
- Pseudotime does not demonstrate developmental ancestry.

## Conclusion

The CD8 result is biologically coherent. It begins from a well-supported naïve/memory-like root and separates into a clear proliferating branch and two related cytotoxic/dysfunction-associated branches. The proliferating lineage is the most distinct; separation between lineages 1 and 3 remains more tentative.
