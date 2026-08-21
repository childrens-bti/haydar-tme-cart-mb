---
title: "Kdm6b-associated T-cell Remodeling: Integrated Summary"
author: "Bicna Song"
date: "2026-08-20"
---

## Scope

This report summarizes the Kdm6b-associated T-cell analysis. The analysis asks
whether T-cell subclusters differ across treatment conditions and whether the
predefined Kdm6b-high state is associated with distinct transcriptional and
immune-program profiles.

## Integrated findings

### 1. Which T-cell subclusters are Kdm6b-high and Kdm6b-low?

The Kdm6b-high group contains:

- C3: NK-like cytotoxic lymphocytes
- C9: γδ T cells
- C12: Th2-like inflammatory T cells

The Kdm6b-low group contains:

- C10: IFN-responsive T cells
- C4: Cycling / proliferating T cells
- C2: Activated effector CD8 T cells (stress)

The remaining seven subclusters are intermediate. These are subcluster-level
definitions based on sample-level Kdm6b expression and are not
individual-cell Kdm6b classifications.

Outputs:

- `results/tcell/tcell_kdm6b_subcluster_ranking_and_groups.tsv`
- `plots/tcell/tcell_kdm6b_subcluster_ranking.pdf`
- `plots/tcell/tcell_kdm6b_violin_by_subcluster.pdf`

### 2. Are Kdm6b-high or Kdm6b-low T-cell subclusters preferentially abundant across conditions?

The abundance pattern is subcluster- and treatment-specific rather than a
uniform enrichment of either Kdm6b group. Kdm6b-high subclusters show selected
enrichment in 41BB-L, including C3 relative to STOP and C12 relative to B7H3
and CD8-41BB. C9 is depleted in 41BB-L relative to CD8-41BB. Kdm6b-low
subclusters also show mixed behavior: C10 and C4 are depleted in 41BB-L in
some comparisons, while C10 is enriched relative to CD8-41BB and C4 is
enriched relative to tumor.

These results support subcluster-specific compositional remodeling rather than
a single group-level Kdm6b abundance effect.

Outputs:

- `plots/tcell/tcell_kdm6b_high_low_subcluster_enrichment_across_conditions.pdf`
- `results/tcell/tcell_kdm6b_high_low_subcluster_differential_abundance.tsv`

### 3. Do Kdm6b-high and Kdm6b-low T cells have different transcriptional profiles?

Yes. The paired sample-level pseudobulk DESeq2 comparison identifies strong
transcriptional differences between the predefined groups. Kdm6b-high cells
show higher expression of genes including `Ifngr1`, `Fosl2`, `Klrb1b`, and
`Vps37b`, while Kdm6b-low cells show higher expression of genes including
`Ly6a`, `Coro1a`, and proliferation-associated genes.

Hallmark GSEA shows Kdm6b-high enrichment for TNFA signaling via NF-kB,
inflammatory response, TGF-beta signaling, IL6-JAK-STAT3 signaling, and
related activation programs. Kdm6b-low is strongly enriched for E2F targets,
G2M checkpoint, MYC targets, oxidative phosphorylation, interferon-alpha and
interferon-gamma response, and other cell-cycle or metabolic programs.

Outputs:

- `plots/tcell/tcell_kdm6b_high_vs_low_pseudobulk_gene_heatmap.pdf`
- `plots/tcell/tcell_kdm6b_high_vs_low_hallmark_gsea_dotplot.pdf`
- `results/tcell/tcell_kdm6b_high_vs_low_paired_deseq2.tsv`
- `results/tcell/tcell_kdm6b_high_vs_low_hallmark_gsea.tsv`
- `results/tcell/tcell_kdm6b_high_vs_low_heatmap_genes.tsv`

### 4. Which T-cell immune programs differ between the Kdm6b groups?

The T-cell UCell analysis uses six predefined programs: cytotoxicity,
expansion/proliferation, memory/stem-like, exhaustion/dysfunction, IFN
response, and chemokine responsiveness.

At the sample level, the cytotoxicity score is higher in Kdm6b-high and remains
significant after adjustment. The IFN-response score is higher in Kdm6b-low
and also remains significant after adjustment. Expansion/proliferation has
zero scores in both groups, while memory/stem-like, exhaustion/dysfunction,
and chemokine responsiveness do not show significant adjusted differences.

Individual gene comparisons are consistent with this mixed pattern. `Gzmb`
and `Ccl3` are higher in Kdm6b-high, whereas `Cxcr3`, `Ifit3`, `Ifng`, `Isg15`,
`Lag3`, `Mki67`, and `Havcr2` are higher in Kdm6b-low after adjustment.

The results support selective cytotoxic and interferon differences rather than
a uniform increase across all T-cell immune programs.

Outputs:

- `plots/tcell/tcell_kdm6b_immune_program_gene_high_vs_low_boxplots.pdf`
- `plots/tcell/tcell_kdm6b_immune_program_gene_violin_by_subcluster.pdf`
- `plots/tcell/tcell_kdm6b_immune_program_ucell_high_vs_low_boxplots.pdf`
- `plots/tcell/tcell_kdm6b_immune_program_ucell_violin_by_subcluster.pdf`
- `results/tcell/tcell_kdm6b_immune_program_gene_high_vs_low.tsv`
- `results/tcell/tcell_kdm6b_immune_program_genes_by_subcluster.tsv`
- `results/tcell/tcell_kdm6b_immune_program_ucell_by_sample_group.tsv`
- `results/tcell/tcell_kdm6b_immune_program_ucell_by_subcluster.tsv`
- `results/tcell/tcell_kdm6b_immune_program_ucell_high_vs_low.tsv`

### 5. Overall interpretation

The T-cell results support an associative Kdm6b-linked remodeling pattern. The
Kdm6b-high group is defined by cytotoxic and inflammatory T-cell states and
shows stronger cytotoxicity, whereas the Kdm6b-low group includes IFN-responsive
and cycling T-cell states and shows stronger interferon and cell-cycle
programs. Differential abundance is mixed across individual subclusters, so
the results do not support a uniform treatment-associated expansion of the
Kdm6b-high or Kdm6b-low group.
