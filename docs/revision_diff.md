# Revision diff — original study re-run through the hardened pipeline

**Question:** after the peer-review hardening (nf-test, tuple-based labeling,
BH correction, correlation-aware gene-set test, MultiQC/UMI/versions, etc.),
does re-running the *original* 12-sample yeast dataset through `nf-g4rnaseq`
change any result?

**How this was produced**

```bash
nextflow run . -profile test_full --outdir results_full_revised --skip_report
```

`-profile test_full` feeds the **identical** count matrix
(`assets/test_data/counts.tsv`, 7154 genes × 12 samples) and full S288C
reference that the original bespoke `analysis/` scripts used, so any difference
is attributable to a *method* change, not to different inputs. Run completed
`[SUCCESS] completed=9 failed=0`. (`--skip_report` only skips the Quarto
render, which is a known macOS-arm64 `deno` limitation, not a pipeline step;
all analysis outputs are produced before it.)

> Scope note: this is a **tool-validation** exercise, not a re-interpretation of
> the pilot data. The pilot's known limitations (n=2, mating-type confound,
> single acute timepoint) are properties of *that dataset* and are addressed by
> the future n=3 matched-genotype experiment — not by the tool.

---

## 1. Predictive / factual outputs — IDENTICAL

The refactor of the bespoke scripts into the packaged tool did **not** perturb
any deterministic result.

| Output | Original (`analysis/`) | Revised (`nf-g4rnaseq`) | Δ |
|---|---|---|---|
| DEGs — PIF1Δ vs WT | 141 (42↑ / 99↓) | 141 | identical |
| DEGs — RecQΔ vs WT | 660 (238↑ / 422↓) | 660 | identical |
| DEGs — PIF1Δ vs RecQΔ | 43 (11↑ / 32↓) | 43 | identical |
| DEGs — PDS (all genotypes) & interactions | 0 | 0 | identical |
| Predicted G4 sites | 7614 | 7614 | identical |
| Hybrid-capable genes (PHQS) | 3138 | 3138 | identical |
| Hybrid-specific genes | 3105 | 3105 | identical |
| GO, RecQΔ-down top terms | rRNA processing / ribosome biogenesis | same | identical |

DESeq2 contrasts, pqsfinder G4 calls, the PHQS hybrid model, and GO enrichment
all reproduce exactly. **No biological/predictive call changed.**

---

## 2. Methods that were *intentionally* changed — now more conservative

Only the two findings the peer review flagged as statistically over-sold move,
and they move in the expected direction: **direction preserved, inflated
significance removed.**

### 2a. G4↔DEG enrichment — now BH-corrected (review item A2)

`bin/g4_enrichment.R` now reports `padj` (Benjamini–Hochberg across the
contrast × set × direction family) and restricts the universe to genes that
were actually scored for G4 features (660 → 642 testable DEGs for RecQΔ).

| Contrast / set | Original | Revised |
|---|---|---|
| RecQΔ, high-conf G4, all DEGs | OR 1.34, **raw p 0.033**, *no correction* | OR 1.44, raw p 0.008, **padj 0.078** |
| RecQΔ, high-conf G4, down DEGs | OR 1.37, raw p 0.058 | OR 1.45, raw p 0.026, **padj 0.136** |

→ The RecQΔ high-confidence-G4 enrichment survives as a **suggestive trend**
(padj ≈ 0.08) rather than a "significant" raw-p hit. Honest, not erased.

### 2b. G4 gene-set GSEA — fgsea → `limma::cameraPR` (review items A1/A3)

The old custom/fgsea GSEA permutes genes and ignores inter-gene correlation, so
~3000-gene G4 sets produced size-inflated p-values (down to 2.6e-8). The
replacement is a **competitive, correlation-aware** test (`cameraPR`,
`use.ranks=TRUE`, `inter.gene.cor=0.01`) reporting Direction + median in-set
statistic + BH-FDR.

**Headline contrast — PDS within RecQΔ** (the "PDS represses G4 genes" claim):

| G4 set | fgsea NES, padj | cameraPR Direction, padj |
|---|---|---|
| coding-strand G4 | −1.34, **2.6e-8** | Down, **0.267** |
| proximal high-conf G4 | −1.31, 0.0032 | Down, 0.267 |
| hybrid-specific | −1.21, 8.2e-5 | Down, 0.267 |
| hybrid-capable | −1.21, 8.8e-5 | Down, 0.267 |

Every set keeps its **Down** direction, but none survive correction — the
1e-8…1e-4 p-values were artifacts of set size + gene-level permutation.

**Direction is preserved across every overlapping contrast** (sanity check that
cameraPR isn't simply nulling everything):

| Contrast | fgsea lead set (dir) | cameraPR same set (dir) | agree? |
|---|---|---|---|
| PDS_in_RecQd | coding-strand-G4 (Down) | Down | ✔ |
| PDS_in_WT | hybrid-capable (Down) | Down | ✔ |
| PDS_in_PIF1d | coding-strand-G4 (Down) | Down | ✔ |
| interaction_RecQd | proximal-highconf (Down) | Down | ✔ |

(One near-zero set in PDS_in_WT, proximal-highconf median −0.056, flips its
label — expected noise at effectively zero effect; padj 0.74 either way.)

---

## 3. Bottom line

- **Reproducibility:** the packaged tool reproduces the original study's
  DE, G4, hybrid-G4, and GO results bit-for-bit on the same inputs.
- **Statistical honesty:** the only numbers that changed are the two the review
  flagged. After correlation-aware testing and multiple-testing correction, the
  G4/hybrid gene-set "enrichments" are revealed as **directionally consistent
  but not statistically significant** in this pilot — exactly the kind of
  over-claim a reviewer would catch, now caught by the tool itself.
- **Net effect on the tool:** strictly safer. Nothing that should be stable
  moved; everything that should have been corrected was corrected. The tool is
  ready to run on the future well-powered (n=3, matched-genotype) dataset, where
  a real coordinated G4 effect — if present — would now clear a defensible bar.
