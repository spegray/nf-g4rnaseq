#!/usr/bin/env Rscript
# render_report.R - assemble a concise Markdown summary from pipeline results.
# (A richer parameterized Quarto template can replace this; this base-R version
#  has no extra dependencies so the pipeline always emits a readable report.)
suppressPackageStartupMessages({ library(optparse); library(data.table) })
opt <- parse_args(OptionParser(option_list=list(
  make_option("--deseq2_dir", type="character"),
  make_option("--g4_dir",     type="character", default=NULL),
  make_option("--hybrid_dir", type="character", default=NULL),
  make_option("--gsea_dir",   type="character", default=NULL),
  make_option("--integrate_dir", type="character", default=NULL),
  make_option("--outdir",     type="character", default=".")
)))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
L <- c("# nf-g4rnaseq summary report", "", paste("Generated:", Sys.Date()), "")
tbl <- function(dt) c("", paste(knitr_kable <- capture.output(print(dt)), collapse="\n"), "")
add <- function(...) L <<- c(L, ...)

ds <- file.path(opt$deseq2_dir, "deg_summary.csv")
if (file.exists(ds)) { add("## Differential expression (significant genes per contrast)", "",
  paste(capture.output(print(fread(ds))), collapse="\n"), "") }
if (!is.null(opt$g4_dir)) {
  f <- file.path(opt$g4_dir, "gene_g4_features.tsv")
  if (file.exists(f)) { g <- fread(f); add(sprintf("## G4 landscape\n\n- %.1f%% of genes have a proximal G4; %.1f%% high-confidence.\n",
      100*mean(g$has_g4), 100*mean(g$has_g4_highconf))) }
}
if (!is.null(opt$hybrid_dir)) {
  f <- file.path(opt$hybrid_dir, "gene_hybrid_g4.tsv")
  if (file.exists(f)) { h <- fread(f); add(sprintf("## RNA:DNA hybrid-G4 (PHQS)\n\n- hybrid-capable: %.1f%%; hybrid-specific: %.1f%%.\n",
      100*mean(h$hybrid_capable), 100*mean(h$hybrid_specific))) }
}
if (!is.null(opt$gsea_dir) && dir.exists(opt$gsea_dir)) {
  gf <- list.files(opt$gsea_dir, pattern="^gseaG4_.*csv$", full.names=TRUE)
  if (length(gf)) { add("## G4 gene-set GSEA (NES; <0 = down-regulated by the contrast)", "")
    for (x in gf) { d <- fread(x); add(paste0("**", sub("gseaG4_|\\.csv","",basename(x)), "**: ",
        paste(sprintf("%s NES=%.2f (padj %.2g)", d$ID, d$NES, d$p.adjust), collapse="; ")), "") } }
}
if (!is.null(opt$integrate_dir)) {
  f <- file.path(opt$integrate_dir, "master_gene_table.tsv")
  if (file.exists(f)) add(sprintf("## Integrated table\n\n- master_gene_table.tsv: %d genes (DE + G4 + hybrid features).\n", nrow(fread(f))))
}
writeLines(L, file.path(opt$outdir, "summary_report.md"))
cat(">>> summary_report.md written\n")
