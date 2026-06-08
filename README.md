# nf-g4rnaseq

**A configurable Nextflow pipeline for RNA-seq differential expression integrated with G-quadruplex (G4) and co-transcriptional RNA:DNA hybrid-G4 (PHQS) analysis.**

[![CI](https://github.com/OWNER/nf-g4rnaseq/actions/workflows/ci.yml/badge.svg)](https://github.com/OWNER/nf-g4rnaseq/actions/workflows/ci.yml)
[![Docker test](https://github.com/OWNER/nf-g4rnaseq/actions/workflows/docker-test.yml/badge.svg)](https://github.com/OWNER/nf-g4rnaseq/actions/workflows/docker-test.yml)
[![Container](https://github.com/OWNER/nf-g4rnaseq/actions/workflows/build-container.yml/badge.svg)](https://github.com/OWNER/nf-g4rnaseq/actions/workflows/build-container.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04-23aa62.svg)](https://www.nextflow.io/)

> Replace `OWNER` in the badge links with your GitHub username/org after publishing. The **Docker test** badge stays green and skips until you publish the image (push a `v*` tag); after that it runs `-profile test,docker` end-to-end on each image build.

Built to answer questions of the form: *does a treatment and/or a genetic background change transcription, and is that change linked to G-quadruplexes (canonical DNA G4s and/or co-transcriptional hybrid G4s)?* You supply a **sample sheet** (samples, replicates, and any experimental factors such as genotype/treatment), a **contrasts sheet**, and a reference; the pipeline produces differential expression, a genome-wide G4 map, a hybrid-G4 (PHQS) map, G4↔DE enrichment, GO/KEGG, threshold-free GSEA (including the genotype×treatment interaction), an integrated per-gene table, and a summary report.

> Origin: generalised from a study of Pyridostatin and the RecQ/PIF1-family G4-resolving helicases in *S. cerevisiae*. The defaults reproduce that study; every organism/design-specific choice is now a parameter.

## Quick start

```bash
# 1. install Nextflow (>=23.04) and a container engine (Docker/Singularity) or conda
# 2. prepare a sample sheet and a contrasts sheet (see assets/ for examples)
nextflow run OWNER/nf-g4rnaseq -r main -profile docker \
    --input        samplesheet.csv \
    --contrasts    contrasts.csv \
    --fasta        genome.fa \
    --gff          annotation.gff3 \
    --orgdb        org.Sc.sgd.db \
    --gene_id_keytype ORF \
    --reference_levels 'genotype=WT;treatment=untreated' \
    --outdir       results
```

Self-contained smoke test (downsampled chromosome I, ~1–2 min — bundled data, nothing else needed):
```bash
nextflow run . -profile test,docker --outdir results
```
Reproduce the full original yeast study (needs the full reference + count matrix):
```bash
nextflow run . -profile test_full,docker --outdir results
```

## Three entry points (auto-detected from the sample sheet, or set `--step`)

| `--step`  | Start from | Sample-sheet columns |
|-----------|-----------|----------------------|
| `fastq`   | raw reads | `sample,fastq_1[,fastq_2],<factors...>` → fastp → STAR → featureCounts |
| `bam`     | aligned (ideally dedup) BAMs | `sample,bam,<factors...>` → featureCounts |
| `counts`  | a gene-count matrix (`--count_matrix`) | `sample,<factors...>` |

## What you get (`results/`)

`reference/` harmonised genome + annotation + SAF · `counts/` count matrix · `deseq2/` per-contrast results + QC (PCA, dispersion, volcano/MA) · `g4/` G4 sites (BED/TSV) + per-gene features + DE↔G4 enrichment · `g4/gene_hybrid_g4.tsv` hybrid-G4 map · `go/` GO+KEGG · `gsea/` GSEA tables + figures · `report/master_gene_table.tsv` + `summary_report.md`.

## Key inputs

- **Sample sheet** (`--input`): one row per sample; a `sample` column plus your factor columns. The `bam`/`fastq_1` paths may be absolute or relative to the launch dir. See `assets/samplesheet_bam_example.csv` (placeholder paths). The bundled `test_bam`/`test_full` regression profiles read `assets/samplesheet_bam_local.csv` — a gitignored copy with your real local paths; create it by copying the example and filling in the `bam` column.
- **Contrasts** (`--contrasts`): `id,type,factor,target,reference,within` where `type` ∈ {`simple`, `pairwise`, `interaction`} — covers within-factor effects, between-level contrasts, and genotype×treatment interactions without editing code. See `assets/contrasts_example.csv`.
- **Reference**: `--fasta` + `--gff`; seqnames are auto-reconciled (`--harmonize length|map|none`).
- **Organism** (optional, enables GO/GSEA): `--orgdb` (Bioconductor OrgDb) + `--gene_id_keytype` (+ `--kegg_organism`).
- **Flag genes** (optional): `--flag_genes` CSV (`gene_id,category`) to mark/exclude confounders (e.g. mating-type or strain-construction genes).

Full parameter list: [`docs/parameters.md`](docs/parameters.md) · usage details: [`docs/usage.md`](docs/usage.md) · outputs: [`docs/output.md`](docs/output.md).

## Reproducibility & the container

`-profile docker` / `-profile singularity` (standard tools via biocontainers; the R analysis via one image) or `-profile conda` (`containers/environment.yml`, no image needed — also what CI uses).

The analysis image is **built and pushed automatically by GitHub Actions** when you push a version tag — no local Docker required:
```bash
git tag v0.1.0 && git push origin v0.1.0      # -> ghcr.io/<owner>/nf-g4rnaseq-r:0.1.0
```
Then make the GHCR package public and set `params.r_container` (or pass `--r_container`). The environment is verified to solve on linux-64, so the build is reliable. See [`containers/README.md`](containers/README.md).

## Citing

If you use nf-g4rnaseq, please cite it — GitHub shows a **"Cite this repository"** button generated from [`CITATION.cff`](CITATION.cff) — and the tools it wraps ([`CITATIONS.md`](CITATIONS.md)), notably DESeq2, pqsfinder, clusterProfiler/fgsea, and, for the RNA:DNA hybrid-G4 (PHQS) model, Zheng et al. 2013, *Nucleic Acids Research* ([10.1093/nar/gkt264](https://doi.org/10.1093/nar/gkt264)).

**Get a citable DOI:** connect the repo to [Zenodo](https://zenodo.org) and publish a GitHub release — Zenodo mints a DOI (metadata pre-filled by [`.zenodo.json`](.zenodo.json)). Then add that DOI to `CITATION.cff`, `CHANGELOG.md`, `nextflow.config` (`manifest.doi`), and as a badge here:
`[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)`

Licensed under [MIT](LICENSE). Changes are tracked in [`CHANGELOG.md`](CHANGELOG.md).
