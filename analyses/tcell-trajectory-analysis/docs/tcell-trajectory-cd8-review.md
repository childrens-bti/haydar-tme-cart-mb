# CD8-like shared trajectory review

Review date: 2026-09-29

Consistency check: 2026-10-01

Label revision: 2026-10-05

Scope: descriptive review of the existing shared `all_conditions` CD8-like results from `03-tcell-cd4-cd8-slingshot.Rmd`. No analysis was rerun. The labels for lineages 1 and 3 were revised on 2026-10-05 using the existing stage 04 tradeSeq lineage rankings (`results/cd4-cd8-tradeseq/tcell_cd8_tradeseq_lineage_rankings.tsv`); fitted changes below are mean condition-specific log2 fitted end/start values.

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
| 1: cytotoxic effector trajectory | Cluster 6 to cluster 1 | Largest tradeSeq increases are cytotoxic genes: `Gzmc` (+7.8 log2), `Gzmb` (+6.8), `Prf1` (+4.1), and `Klrg1` (+3.4); no significant change in `Tox` (FDR 0.28); cytotoxic AUCell activity stays high at the endpoint; endpoint is 59.0% ProjecTILs `CD8_Tex`; median pseudotime increases from 9.4 to 33.2 | Inhibitory receptors (`Pdcd1`, `Havcr2`, `Ctla4`) also rise, and ProjecTILs labels the endpoint mostly `CD8_Tex`; the endpoint also has NKT/NK-like SingleR support |
| 2: proliferating effector/dysfunction trajectory | Cluster 6 to cluster 5 | `Mki67`, `Top2a`, proliferation AUCell, cytotoxicity, and exhaustion-associated programs increase; median pseudotime increases from 8.8 to 46.6 | Cell cycle may drive much of this branch's separation |
| 3: exhausted/dysfunctional trajectory | Cluster 6 to cluster 3 | Largest tradeSeq increases are exhaustion-associated genes: `Xcl1` (+4.4 log2), `Pdcd1` (+4.2), `Lag3` (+3.2), and `Tox` (+1.7, FDR 0.06); no significant change in `Prf1` (FDR 1.0); `Gzmb` and `Ccl5` peak earlier in the stage 03 curves; exhaustion AUCell keeps rising to the endpoint; endpoint is 43.8% ProjecTILs `CD8_Tex`; median pseudotime increases from 8.8 to 27.6 | Shares much of its early route with lineage 1; `Gzmb` still increases (+2.3 log2); the endpoint is mixed and not a pure exhausted population |

## Condition representation

All six conditions have cells with non-missing pseudotime on all three lineages, but their representation is unequal. B7H3 contributes many cells with lineage 1 pseudotime, while CD8-41BB contributes many cells with lineage 3 pseudotime. Cells can contribute to multiple lineages through Slingshot weights, so these are not exclusive fate counts. The condition-specific smooth curves are descriptive, not evidence of condition-specific fate preference.

## Main limitations

- The proliferating branch may be partly defined by cell cycle.
- Lineage 1 may contain an NK/NKT-like cytotoxic population.
- The exhaustion AUCell score declines near the endpoint of lineages 1 and 2 in several constructs; this is not explained by the current analysis, but is consistent with lineage 1 retaining a cytotoxic program.
- Lineages 1 and 3 share much of their early route; they differ mainly in which program dominates at the endpoint (cytotoxic for lineage 1, exhaustion-associated for lineage 3).
- The stage 03 gene curves are descriptive smoothers without uncertainty intervals or significance tests; the stage 04 tradeSeq fitted changes and association tests are cell-level and exploratory.
- Unequal cell recovery and only two samples per condition limit condition-level interpretation.
- Pseudotime does not demonstrate developmental ancestry.

## Conclusion

The CD8 result is biologically coherent. It begins from a well-supported naïve/memory-like root and separates into a clear proliferating branch, a cytotoxic effector branch (lineage 1), and an exhausted/dysfunctional branch (lineage 3). The proliferating lineage is the most distinct. Lineages 1 and 3 share much of their early route and both carry some cytotoxic and exhaustion-associated features, so their separation is based on which program dominates at the endpoint and remains more tentative.
