#!/usr/bin/env Rscript
# =============================================================================
# integrate.R  -  Merge DE + G4 + hybrid-G4 features into one master table
# Generalises analysis/scripts/08_integrate.R.
# =============================================================================
suppressPackageStartupMessages({ library(optparse); library(data.table) })
opt <- parse_args(OptionParser(option_list=list(
  make_option("--gene_info",   type="character"),
  make_option("--g4_features", type="character", default=NULL),
  make_option("--hybrid",      type="character", default=NULL),
  make_option("--deseq2_dir",  type="character", help="dir containing res_list.rds"),
  make_option("--flag_genes",  type="character", default=NULL),
  make_option("--padj",        type="double", default=0.05),
  make_option("--outdir",      type="character", default=".")
)))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
gi <- fread(opt$gene_info)
rl <- readRDS(file.path(opt$deseq2_dir, "res_list.rds"))
M <- gi[, .(gene_id, std_name, orf_class, biotype, chr, start, end, strand)]
# baseMean from the first contrast
b <- as.data.table(rl[[1]])[, .(gene_id, baseMean=round(baseMean,1))]
M <- merge(M, b, by="gene_id", all.x=TRUE)
# flags
if (!is.null(opt$flag_genes) && file.exists(opt$flag_genes)) {
  fl <- fread(opt$flag_genes); for (cat in unique(fl$category)) M[[paste0("is_",cat)]] <- M$gene_id %in% fl[category==cat]$gene_id
}
# per-contrast LFC + padj
for (cn in names(rl)) {
  d <- as.data.table(rl[[cn]])[, .(gene_id, lfc=round(log2FoldChange,3), padj=signif(padj,3))]
  setnames(d, c("lfc","padj"), paste0(c("LFC_","padj_"), cn)); M <- merge(M, d, by="gene_id", all.x=TRUE)
}
if (!is.null(opt$g4_features) && file.exists(opt$g4_features)) M <- merge(M, fread(opt$g4_features), by="gene_id", all.x=TRUE)
if (!is.null(opt$hybrid) && file.exists(opt$hybrid))
  M <- merge(M, fread(opt$hybrid)[, .(gene_id, sense_G_frac, sense_G3_tracts, hybrid_capable, hybrid_specific, expressed)], by="gene_id", all.x=TRUE)
setcolorder(M, c("gene_id","std_name","biotype","chr","start","end","strand","baseMean"))
fwrite(M[order(chr,start)], file.path(opt$outdir,"master_gene_table.tsv"), sep="\t")
cat(sprintf(">>> master_gene_table.tsv: %d genes x %d columns\n", nrow(M), ncol(M)))
cat(">>> integrate complete.\n")
