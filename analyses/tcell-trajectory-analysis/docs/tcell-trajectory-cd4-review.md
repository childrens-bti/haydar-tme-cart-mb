# CD4-like shared trajectory review

Review date: 2026-09-29

Consistency check: 2026-10-01

Scope: descriptive review of the existing shared `all_conditions` CD4-like results from `03-tcell-cd4-cd8-slingshot.Rmd`. No analysis was rerun.

## Purpose

The trajectory orders CD4-like cells along connected transcriptional states to identify possible regulatory, proliferating, and activated helper-like branches. Pseudotime is an inferred state ordering, not direct observation of cells changing over time.

## Overview

- 2,505 cells, six retained clusters, and three lineages.
- Script 02 contains no CD4 cluster annotated as naïve/memory by ProjecTILs.
- Cluster 3 was selected as the root using the highest available naïve/memory marker score.
- Cluster 3 is 62.3% ProjecTILs `Th1` and is most often previously annotated as quiescent T cells.

The root selector first searches for a ProjecTILs label containing `naive`. Because the CD4 results contain only `Th1` and `Treg` cluster labels, it used its fallback: select the cluster with the highest positive median score across `Tcf7`, `Lef1`, `Sell`, `Ccr7`, and `Il7r`. Cluster 3 ranked first with a median score of 0.492; the next-highest cluster was cluster 5 at 0.348. Cluster 3 was also the only cluster whose dominant previous annotation was quiescent T cells.

Cluster 3 is therefore the **least-differentiated available CD4 state** and is appropriately described as a quiescent Th1-like origin. Its selection does not imply that it is a canonical naïve CD4 population; no such population is present in these data.

## Lineages

| Lineage | Main result | Supporting evidence | Main caution |
|---|---|---|---|
| 1: regulatory/Treg-associated | Cluster 3 to cluster 2 | Cluster 2 is 81.9% ProjecTILs `Treg`, 91.9% previously annotated Treg, and Treg-like by SingleR; median pseudotime increases from 20.4 to 38.7 | Core markers such as `Foxp3`, `Il2ra`, `Ikzf2`, and `Tigit` were not plotted, so the endpoint is better supported than the full path |
| 2: proliferating CD4 | Cluster 3 to cluster 1 | Cluster 1 is 75.9% cycling/proliferating; `Mki67`, `Top2a`, and proliferation AUCell activity rise strongly; median pseudotime increases from 19.4 to 53.0 | Cell cycle may drive much of the branch geometry; the helper identity of the cycling cells is unresolved |
| 3: activated helper-like | Cluster 3 to cluster 5 | Cluster 5 is 81.9% ProjecTILs `Th1` and 90.6% previously activated helper-like; `Ifng` and `Pdcd1` rise late; median pseudotime increases from 20.5 to 41.2 | Current markers do not clearly distinguish productive helper activation from chronic activation or dysfunction |

## Condition representation

All six conditions have cells with non-missing pseudotime on all three lineages, but their representation is unequal. CD8-41BB contributes the most cells to each lineage, while tumor contributes only 26 cells to lineage 2. Cells can contribute to multiple lineages through Slingshot weights, so these are not exclusive fate counts. The condition-specific smooth curves are descriptive, not evidence of condition differences.

## Main limitations

- No canonical naïve/memory CD4 cluster is present, so progression is relative to the least-differentiated available quiescent Th1-like state.
- The plotted gene panel lacks important CD4-specific regulatory and helper markers.
- The proliferating branch may be dominated by cell cycle.
- Unequal cell recovery and only two samples per condition limit condition-level interpretation.
- Pseudotime does not demonstrate developmental ancestry.

## Conclusion

The CD4 result begins from the best-supported quiescent Th1-like state available in this compartment and separates toward a strongly supported Treg-associated cluster, a clear proliferating branch, and an activated Th1/helper-like branch. It should be interpreted as progression relative to that observed origin, not as a naïve-to-effector trajectory.
