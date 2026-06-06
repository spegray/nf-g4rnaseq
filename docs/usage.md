# nf-g4rnaseq — usage

## 1. Sample sheet (`--input`)
One row per sample. A `sample` column (unique IDs) plus any experimental-factor
columns referenced in `--design_formula` (e.g. `genotype`, `treatment`,
`replicate`). Columns depend on the entry point:

| Entry | Required columns |
|-------|------------------|
| FASTQ | `sample,fastq_1[,fastq_2],<factors...>` |
| BAM   | `sample,bam,<factors...>` |
| counts| `sample,<factors...>` (+ `--count_matrix file.tsv`) |

Extra factor columns are allowed and may be used in the design/contrasts.
The entry point is auto-detected from the columns (`fastq_1` → fastq, `bam` →
bam, else counts) or set explicitly with `--step`.

## 2. Contrasts (`--contrasts`)
CSV `id,type,factor,target,reference,within`. One row per comparison:

| `type` | Meaning | Uses columns |
|--------|---------|--------------|
| `simple` | effect of `factor` (target vs reference) **within** a level of another factor | `factor,target,reference,within` (e.g. `within=genotype=WT`) |
| `pairwise` | `factor` target vs reference (main effect) | `factor,target,reference` |
| `interaction` | does the `f1` effect differ across `f2`? | `factor=f1:f2,target=l1:l2,reference=lr1:lr2` |

Example (the originating study):
```csv
id,type,factor,target,reference,within
PDS_in_WT,simple,treatment,PDS,untreated,genotype=WT
RecQd_vs_WT,pairwise,genotype,RecQd,WT,
interaction_RecQd,interaction,genotype:treatment,RecQd:PDS,WT:untreated,
```
Set `--design_formula '~ replicate + genotype * treatment'` and
`--reference_levels 'genotype=WT;treatment=untreated'` so the reference levels
match the contrasts. `simple` effects use apeglm shrinkage (relevel+refit);
`pairwise` use ashr; `interaction` uses the interaction coefficient.

## 3. Reference & seqname harmonization
`--fasta genome.fa --gff annotation.gff3`. The genome, annotation, and (BAM
entry) the reads may use different chromosome names; `--harmonize` reconciles them:
- `length` (default): match contigs by length to the read coordinate system.
- `map`: use an explicit `--seqname_map old<TAB>new` table.
- `none`: assume they already match.

## 4. Organism (enables GO/KEGG/GSEA-GO)
`--orgdb org.Sc.sgd.db` (any Bioconductor OrgDb), `--gene_id_keytype ORF`
(the keytype matching your gene IDs, e.g. `ENSEMBL`), optional
`--kegg_organism sce` (KEGG code; needs internet). If `--orgdb` is omitted, the
DE + G4 + hybrid core still runs (GO/GSEA-GO are skipped); custom G4 gene-set
GSEA still runs.

## 5. Optional confounder flagging
`--flag_genes flags.csv` (`gene_id,category`) marks genes (e.g. mating-type,
strain-construction) so they are excluded from GO/enrichment and annotated in
the master table.

## 6. Profiles & running
```bash
# containerized (recommended for others / HPC / cloud)
nextflow run OWNER/nf-g4rnaseq -r main -profile docker --input ... --outdir results
# or singularity / conda
-profile singularity      -profile conda
# reproduce the bundled yeast example
nextflow run . -profile test,docker
```
Resume after an interruption with `-resume`. Tune resources via `--max_cpus`,
`--max_memory`. Skip stages with `--skip_g4 --skip_hybrid --skip_go --skip_gsea --skip_report`.

## UMI libraries (FASTQ entry)
The FASTQ path runs fastp → STAR → featureCounts. For UMI-based libraries,
insert `umi_tools extract` (before STAR) and `umi_tools dedup` (after); the
`umi_tools` tool is already in the environment. See `modules/local/align.nf`.
