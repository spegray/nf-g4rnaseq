# nf-g4rnaseq — outputs (`results/`)

| Path | What it is |
|------|-----------|
| `reference/` | genome + annotation with reconciled chromosome names; `gene_info.tsv` (gene_id, std_name, coords, biotype); `genes.saf` (featureCounts annotation) |
| `counts/counts.tsv` | gene × sample raw count matrix (FASTQ/BAM entry); `*.summary` = per-sample assignment QC |
| `deseq2/` | per-contrast results `*.csv` (shrunken log2FC, Wald p/padj, gene annotation, flags); `deg_summary.csv` (sig genes per contrast); `dds.rds`, `res_list.rds`; `figures/` = PCA, dispersion, sample-distance heatmap, MA + volcano per contrast |
| `g4/g4_sites.{tsv,bed}` | every predicted G4: position, strand, pqsfinder score, #tetrads (≥3 = canonical, 2 = two-quartet), bulges/mismatches. BED opens in IGV |
| `g4/gene_g4_features.tsv` | per-gene proximal-G4 counts (promoter+body), split by coding vs template strand and by stringency (all / canonical / high-confidence) |
| `g4/gene_hybrid_g4.tsv` | per-gene RNA:DNA hybrid-G4 (PHQS) scan: coding-strand G-tract clusters; `hybrid_capable` (≥2 tracts) and `hybrid_specific` (2–3 tracts → forms a G4 only co-transcriptionally) |
| `g4/g4_deg_enrichment.csv` | Fisher tests: are DEGs (all/up/down, per contrast) enriched for proximal/hybrid G4s? (with an optional ribosome-excluded robustness set) |
| `go/` | GO (BP/MF/CC) and KEGG over-representation per contrast (up & down), with dotplots — only if `--orgdb` given |
| `gsea/` | threshold-free GSEA per contrast: `gseGO_*.csv` (GO sets) and `gseaG4_*.csv` (custom G4 gene sets); `figures/G4_GSEA_NES_byContrast.png`. Works on simple, pairwise AND interaction contrasts (NES<0 = the set is coordinately down-regulated) |
| `integrate/master_gene_table.tsv` | one row per gene: expression, per-contrast log2FC+padj, proximal-G4 features, hybrid-G4 class, flags — the integrated backbone |
| `report/report.html` | the full **Quarto** report (narrative + tables + embedded figures); rendered where Quarto is available (bundled in the container) |
| `report/summary_report.md` | dependency-free Markdown summary, always produced (and the fallback if Quarto is unavailable) |

## Interpreting the G4 columns
- **canonical G4** = ≥4 G-tracts on one DNA strand (`n_tetrads`≥3, classic motif).
- **hybrid-capable / hybrid-specific** = the co-transcriptional model (Zheng et al. 2013): the nascent transcript donates G-tracts so a gene with only 2–3 coding-strand tracts — too few for a DNA-only G4 — can still form a G4 during transcription. `hybrid_specific` genes are the novel, transcription-dependent class.
- In GSEA, a **negative NES for a G4 gene set** means that set tends to be down-regulated in that contrast (e.g. a G4-stabilising drug repressing G4-bearing genes).

## Caveats (carried from the study; apply generally)
- Low replication limits per-gene power; GSEA aggregates across genes and is more robust. 
- Permissive G4 calling has low specificity — prefer the high-confidence subset for enrichment claims.
- Confounded design factors (e.g. a factor perfectly aligned with a batch or, in the study, mating type) should be flagged via `--flag_genes` and interpreted with care; prefer contrasts that hold confounders constant.
