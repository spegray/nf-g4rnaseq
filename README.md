# nf-g4rnaseq

**A configurable Nextflow pipeline for RNA-seq differential expression integrated with G-quadruplex (G4) and co-transcriptional RNA:DNA hybrid-G4 (PHQS) analysis.**

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

Reproduce the original yeast study end-to-end:
```bash
nextflow run . -profile test,docker --outdir results   # uses assets/*_example.csv
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

- **Sample sheet** (`--input`): one row per sample; a `sample` column plus your factor columns. See `assets/samplesheet_bam_example.csv`.
- **Contrasts** (`--contrasts`): `id,type,factor,target,reference,within` where `type` ∈ {`simple`, `pairwise`, `interaction`} — covers within-factor effects, between-level contrasts, and genotype×treatment interactions without editing code. See `assets/contrasts_example.csv`.
- **Reference**: `--fasta` + `--gff`; seqnames are auto-reconciled (`--harmonize length|map|none`).
- **Organism** (optional, enables GO/GSEA): `--orgdb` (Bioconductor OrgDb) + `--gene_id_keytype` (+ `--kegg_organism`).
- **Flag genes** (optional): `--flag_genes` CSV (`gene_id,category`) to mark/exclude confounders (e.g. mating-type or strain-construction genes).

Full parameter list: [`docs/parameters.md`](docs/parameters.md) · usage details: [`docs/usage.md`](docs/usage.md) · outputs: [`docs/output.md`](docs/output.md).

## Reproducibility

`-profile docker` / `-profile singularity` (standard tools via biocontainers; the R analysis via the image built from `containers/Dockerfile`) or `-profile conda` (`containers/environment.yml`). Pin versions on release.

## Citation

If you use nf-g4rnaseq, please cite the tools it wraps (see [`CITATIONS.md`](CITATIONS.md)), notably DESeq2, pqsfinder, clusterProfiler/fgsea, and — for the hybrid-G4 (PHQS) model — Zheng et al. 2013, *Nucleic Acids Research*. Licensed under MIT.
