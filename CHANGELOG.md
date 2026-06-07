# Changelog

All notable changes to **nf-g4rnaseq** are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versioning follows [Semantic Versioning](https://semver.org/).

## [Unreleased] — peer-review hardening
### Fixed
- **Sample labeling**: count columns are now labeled via a sorted `[sample,bam]`
  tuple list (not two independent `.collect()`s), removing a silent mislabel risk;
  added a column-count assertion.
- **Gene-set statistics**: replaced size-inflated fgsea on the large G4 gene sets
  with `limma::cameraPR` (competitive, correlation-aware, calibrated for any set
  size); added Benjamini-Hochberg correction across the G4↔DE Fisher family.
- **FASTQ portability**: STAR now reads decompressed FASTQ (the `--readFilesCommand`
  spawn fails on some platforms).
### Added
- MultiQC QC aggregation; `software_versions.yml` provenance; optional UMI handling
  (`--umi`/`--umi_pattern`); `--exclude_dubious` gene-model option; a low-replication
  (<3/group) warning; strandedness-probe transparency; `test_bam`/`test_fastq`
  profiles with a bundled tiny FASTQ dataset.
### Validated
- BAM entry end-to-end on real BAMs (auto strandedness = 1; reproduced DEG counts;
  HXT13 labeling spike correct). Contrast engine confirmed on a 1-factor/3-level
  design (generality). FASTQ wiring validated locally except STAR's macOS build;
  exercised on the Linux CI runner.

## [0.1.0] - 2026-06-06
First release — a configurable Nextflow (DSL2) pipeline generalised from a study
of Pyridostatin and the RecQ/PIF1-family G4-resolving helicases in *S. cerevisiae*.

### Added
- **Entry points**: FASTQ (fastp → STAR → featureCounts), BAM (featureCounts),
  and count-matrix — auto-detected from the sample sheet.
- **Reference harmonization** of genome/annotation/BAM chromosome names
  (`length` / `map` / `none`); builds the featureCounts SAF + gene_info.
- **Differential expression** (DESeq2) with a config-driven design formula and a
  contrasts engine supporting `simple` (within-factor), `pairwise`, and
  `interaction` contrasts; apeglm/ashr shrinkage; PCA/dispersion/MA/volcano QC.
- **G-quadruplex prediction** (pqsfinder, both strands, scored) and per-gene
  proximal-G4 features (promoter + body, strand-aware, stringency tiers).
- **RNA:DNA hybrid-G4 (PHQS) scan** implementing the Zheng et al. 2013 model.
- **G4 ↔ DE enrichment** (directional, with a ribosome-excluded robustness set),
  **GO/KEGG** (clusterProfiler), and **GSEA** including the genotype×treatment
  interaction and custom G4 gene sets.
- **Integrated** per-gene master table and a parameterized **Quarto report**
  (HTML, with a dependency-free Markdown fallback).
- **Organism-configurable** (OrgDb / KEGG / seqname harmonization); optional
  confounder flagging via `--flag_genes`.
- **Reproducibility**: container (Docker/Singularity) + conda profiles; pinned
  `environment.yml` verified to solve on linux-64.
- **Test data & CI**: bundled downsampled (chromosome I) dataset; self-contained
  `-profile test` (~30 s) and full-data `-profile test_full`; CI for compile,
  end-to-end, automated container build/push to GHCR, and a gated
  `-profile test,docker` test.

### Validated
- Reproduces the originating study **exactly**: RecQΔ 660 / PIF1Δ 141 /
  PIF1Δ-vs-RecQΔ 43 DEGs; 7,614 pqsfinder G4 sites; 3,105 hybrid-specific genes;
  PDS-in-RecQΔ coding-strand-G4 GSEA NES −1.34 (padj 2.6e-8).

[Unreleased]: https://github.com/OWNER/nf-g4rnaseq/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/OWNER/nf-g4rnaseq/releases/tag/v0.1.0
