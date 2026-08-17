---
title: "Kdm6b-associated Myeloid Remodeling: Integrated Summary"
author: "Bicna Song"
date: "2026-08-17"
---

## Scope

This report provides an integrated summary of the Kdm6b-associated myeloid
remodeling analysis. The analysis evaluates whether Kdm6b-associated myeloid
subclusters differ across treatment conditions and whether the Kdm6b-high state
is associated with immune-supportive programs.

## Integrated findings

### 1. Is Kdm6b expression highest in the effective CAR condition?

Yes, descriptively. 41BB-L has the highest median sample-level pseudobulk
Kdm6b expression in all three compartments:

- All cells: 15.81
- Myeloid cells: 15.13
- T cells: 10.63

Output:

- `plots/myeloid/kdm6b_pseudobulk_expression_by_compartment.pdf`

### 2. Is the Kdm6b signal primarily associated with myeloid cells?

The data support focusing downstream analyses on the myeloid compartment, but
they do not establish a quantitatively stronger myeloid signal by directly
comparing y-axis values between compartments. Each compartment was
normalized independently, so absolute values across all cells, T cells, and
myeloid cells are not directly comparable.

Within the myeloid compartment, Kdm6b expression is heterogeneous across the
UMAP and subclusters, and the myeloid-specific ranking provides the basis for
the downstream state analysis. Kdm6b is also detectable in T cells, so the
results do not support describing it as myeloid-exclusive.

Outputs:

- `plots/myeloid/kdm6b/kdm6b_myeloid_featureplot.pdf`
- `plots/myeloid/kdm6b/kdm6b_myeloid_violin_by_subcluster.pdf`
- `results/myeloid/kdm6b/kdm6b_pseudobulk_by_sample.tsv`

### 3. Which myeloid subclusters are Kdm6b-high?

The four Kdm6b-high subclusters are:

- C4: Inflammatory activated microglia
- C10: Perivascular macrophages (chemokine-producing)
- C0: Early activated microglia
- C11: Antigen-presenting monocytes

The four Kdm6b-low subclusters are C8, C13, C2, and C9. The remaining eight
subclusters are intermediate. These are subcluster-level definitions based on
median sample-level mean log-normalized Kdm6b expression; they are not
individual-cell Kdm6b classifications.

Output:

- `results/myeloid/kdm6b/kdm6b_subcluster_ranking_and_groups.tsv`
- `plots/myeloid/kdm6b/kdm6b_subcluster_ranking.pdf`

### 4. Are Kdm6b-high subclusters enriched in effective CAR conditions?

Some are, but the pattern is subcluster- and treatment-specific rather than a
uniform expansion of the Kdm6b-high group.

Using the all-16-subcluster sccomp model with 41BB-L as the reference:

- Relative to B7H3, three Kdm6b-high subclusters are enriched in 41BB-L and one
  is enriched in B7H3.
- Relative to CD28-41BB, no Kdm6b-high subcluster is significantly different.
- Relative to CD8-41BB, one Kdm6b-high subcluster is enriched in CD8-41BB.
- Relative to STOP, one Kdm6b-high subcluster is enriched in STOP.
- Relative to tumor, one Kdm6b-high subcluster is enriched in 41BB-L and two are
  enriched in tumor.

The corresponding Kdm6b-low results are also mixed. For example, C13 is
depleted in 41BB-L across comparisons, while C9 is enriched in 41BB-L in
several comparisons. Therefore, the abundance results support
subcluster-specific compositional remodeling rather than a single group-level
Kdm6b-high enrichment effect.

Outputs:

- `plots/myeloid/kdm6b/kdm6b_high_low_subcluster_enrichment_across_conditions.pdf`
- `results/myeloid/kdm6b/kdm6b_high_low_subcluster_differential_abundance.tsv`

### 5. Are Kdm6b-low subclusters depleted or unchanged in 41BB-L?

Neither conclusion applies uniformly. Kdm6b-low subclusters show both
enrichment and depletion depending on the subcluster and comparison. The
results should therefore be interpreted at the individual-subcluster level,
not as evidence that the entire Kdm6b-low compartment is depleted or unchanged
in 41BB-L.

### 6. Do Kdm6b-high subclusters show immune-supportive programs?

The evidence is strongest for inflammatory and interferon-related activity,
with additional evidence for antigen-presentation and selected chemokine genes.
It is not uniform across every gene or scoring method.

#### Transcriptional comparison

The paired sample-level pseudobulk DESeq2 comparison identifies strong
differences between the predefined Kdm6b-high and Kdm6b-low compartments. Genes
higher in Kdm6b-high include inflammatory and activation-associated genes such
as `Csf1`, `Rela`, `Cd83`, and `Icam1`. Hallmark GSEA shows Kdm6b-high
enrichment for:

- TNFA signaling via NF-kB
- Inflammatory response
- Interferon-gamma response
- IL2-STAT5 signaling
- IL6-JAK-STAT3 signaling

Kdm6b-low is enriched for oxidative phosphorylation, fatty-acid metabolism,
and related metabolic programs. These comparisons distinguish the selected
subcluster groups; they do not establish that Kdm6b causes the observed
transcriptional programs.

Outputs:

- `plots/myeloid/kdm6b/myeloid_kdm6b_high_vs_low_pseudobulk_gene_heatmap.pdf`
- `plots/myeloid/kdm6b/myeloid_kdm6b_high_vs_low_hallmark_gsea_dotplot.pdf`
- `results/myeloid/kdm6b/myeloid_kdm6b_high_vs_low_paired_deseq2.tsv`
- `results/myeloid/kdm6b/myeloid_kdm6b_high_vs_low_hallmark_gsea.tsv`
- `results/myeloid/kdm6b/myeloid_kdm6b_high_vs_low_heatmap_genes.tsv`

#### Immune-program genes and UCell scores

At the sample level, several individual genes have higher nominal paired values
in Kdm6b-high,
including `Ccl5`, `Ccr7`, `Cxcl9`, `Cxcl10`, `Cd74`, `H2-Ab1`, `H2-Eb1`,
`Il12b`, `Irf1`, `Tnf`, and `Xcr1` after paired testing. `B2m`, `Ifngr1`,
`Psmb8`, and `Stat1` show the opposite direction in this comparison.

The UCell comparison shows a significant increase for the inflammatory
activation module in Kdm6b-high. The antigen-presentation module is not
significantly different, and the T-cell-recruitment module has zero median
scores in both groups in the sample-level comparison. The gene-level and
UCell-level results therefore support inflammatory activation more clearly
than a broad increase across all immune programs.

Outputs:

- `plots/myeloid/kdm6b/myeloid_kdm6b_immune_program_gene_high_vs_low_boxplots.pdf`
- `plots/myeloid/kdm6b/myeloid_kdm6b_immune_program_gene_violin_by_subcluster.pdf`
- `plots/myeloid/kdm6b/myeloid_kdm6b_immune_program_ucell_high_vs_low_boxplots.pdf`
- `plots/myeloid/kdm6b/myeloid_kdm6b_immune_program_ucell_violin_by_subcluster.pdf`
- `results/myeloid/kdm6b/myeloid_kdm6b_immune_program_gene_high_vs_low.tsv`
- `results/myeloid/kdm6b/myeloid_kdm6b_immune_program_genes_by_subcluster.tsv`
- `results/myeloid/kdm6b/myeloid_kdm6b_immune_program_ucell_by_sample_group.tsv`
- `results/myeloid/kdm6b/myeloid_kdm6b_immune_program_ucell_by_subcluster.tsv`
- `results/myeloid/kdm6b/myeloid_kdm6b_immune_program_ucell_high_vs_low.tsv`

#### Sample-level Kdm6b correlations

Across all myeloid cells, Kdm6b expression is positively correlated with genes
from all three program categories. The strongest association is with `Xcl1`
in the T-cell-recruitment program. Several antigen-presentation genes, including
`Tap1`, `Tapbp`, `H2-Eb1`, `H2-Ab1`, and `Cd74`, show positive correlations, as
do inflammatory genes such as `Irf7`, `Stat1`, `Cxcl10`, and `Tnf`.

After multiple-testing adjustment, `Xcl1` is the clearest statistically
supported correlation. Most other genes have positive nominal associations but
adjusted p-values above 0.05. With only 12 samples, these correlations should
be treated as supportive trends rather than definitive evidence.

Outputs:

- `plots/myeloid/kdm6b/myeloid_kdm6b_correlation_antigen_presentation_scatterplots.pdf`
- `plots/myeloid/kdm6b/myeloid_kdm6b_correlation_antigen_presentation_dotplot.pdf`
- `plots/myeloid/kdm6b/myeloid_kdm6b_correlation_inflammatory_activation_scatterplots.pdf`
- `plots/myeloid/kdm6b/myeloid_kdm6b_correlation_inflammatory_activation_dotplot.pdf`
- `plots/myeloid/kdm6b/myeloid_kdm6b_correlation_t_cell_recruitment_scatterplots.pdf`
- `plots/myeloid/kdm6b/myeloid_kdm6b_correlation_t_cell_recruitment_dotplot.pdf`
- `results/myeloid/kdm6b/myeloid_kdm6b_gene_correlation_summary.tsv`
- `results/myeloid/kdm6b/myeloid_kdm6b_gene_mean_expression_by_sample.tsv`

### 7. Does the Kdm6b-high state appear consistent with immune-supportive remodeling?

The combined results are consistent with a Kdm6b-high myeloid state associated
with inflammatory remodeling that could support CAR T-cell activity. This
interpretation is supported by:

- Highest descriptive Kdm6b expression in 41BB-L within each compartment.
- Kdm6b-high subclusters that include inflammatory activated microglia,
  chemokine-producing perivascular macrophages, and antigen-presenting
  monocytes.
- High-versus-low enrichment for TNFA/NF-kB, inflammatory, interferon-gamma,
  IL2-STAT5, and IL6-JAK-STAT3 pathways.
- Positive sample-level correlations with antigen-presentation,
  inflammatory, and selected recruitment genes.
- Higher Kdm6b-high inflammatory UCell scores.

The conclusion is associative and qualified. Differential abundance is mixed
across individual subclusters, the antigen-presentation UCell score is not
significantly different, the T-cell-recruitment UCell score is uninformative at
the sample-group level, and most correlations do not remain significant after
multiple-testing adjustment. The results support a remodeling hypothesis, not a
claim that Kdm6b directly causes an immune-supportive state or that the state
is uniformly enriched in 41BB-L.
