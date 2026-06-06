#!/usr/bin/env Rscript
# =============================================================================
# run_deseq2.R  -  Configurable differential expression (nf-g4rnaseq)
# -----------------------------------------------------------------------------
# Generalises analysis/scripts/03_deseq2.R: instead of a hardcoded design and
# eight hand-written contrasts, it takes a count matrix, a sample sheet, a design
# FORMULA, and a CONTRASTS table, and produces shrunken results + QC for every
# requested contrast. Three contrast types are supported (see --contrasts):
#   simple      effect of <factor> (target vs reference) WITHIN a level of another
#               factor (column `within`, e.g. genotype=WT)  -> relevel & refit,
#               apeglm-shrunken coefficient   (e.g. "PDS within genotype=WT")
#   pairwise    <factor> target vs reference                 -> ashr-shrunken
#   interaction <f1:f2> target<l1:l2> vs reference<l1:l2>     -> apeglm coefficient
#               (does the f1 effect differ across f2? the "exacerbate/alleviate")
# =============================================================================
suppressPackageStartupMessages({
  library(optparse); library(DESeq2); library(data.table)
  library(ggplot2); library(ggrepel); library(pheatmap); library(RColorBrewer)
})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts",      type="character", help="gene x sample count matrix (TSV; col 1 = gene_id)"),
  make_option("--samplesheet", type="character", help="CSV with column 'sample' + factor columns"),
  make_option("--gene_info",   type="character", default=NULL, help="TSV: gene_id, std_name, chr, start, end, strand, biotype"),
  make_option("--design",      type="character", default="~ replicate + genotype * treatment", help="DESeq2 design formula"),
  make_option("--contrasts",   type="character", help="CSV: id,type,factor,target,reference,within"),
  make_option("--reference_levels", type="character", default="", help="'factor=level;factor2=level2' reference levels"),
  make_option("--flag_genes",  type="character", default=NULL, help="optional CSV: gene_id,category (genes to flag, e.g. mating-type/construction)"),
  make_option("--min_count",   type="integer",   default=10L, help="drop genes with < this many total reads"),
  make_option("--padj",        type="double",    default=0.05, help="adjusted-p significance threshold"),
  make_option("--outdir",      type="character", default=".", help="output directory")
)))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
figdir <- file.path(opt$outdir, "figures"); dir.create(figdir, showWarnings=FALSE)
set.seed(1)

# ---- inputs ----------------------------------------------------------------
cm <- fread(opt$counts); gene_ids <- cm[[1]]
counts <- as.matrix(cm[, -1, with=FALSE]); rownames(counts) <- gene_ids; mode(counts) <- "integer"
ss <- fread(opt$samplesheet)
stopifnot("sample" %in% names(ss))
ss <- ss[match(colnames(counts), ss$sample)]               # align metadata to matrix columns
stopifnot(!any(is.na(ss$sample)), identical(as.character(ss$sample), colnames(counts)))

# factor columns = everything in the design formula
design_fm <- as.formula(opt$design)
fct <- all.vars(design_fm)
refs <- list()
if (nzchar(opt$reference_levels))
  for (kv in strsplit(opt$reference_levels, ";")[[1]]) { p <- strsplit(trimws(kv), "=")[[1]]; refs[[p[1]]] <- p[2] }
cd <- as.data.frame(ss)
for (f in fct) {
  cd[[f]] <- factor(cd[[f]])
  if (!is.null(refs[[f]])) cd[[f]] <- relevel(cd[[f]], ref = refs[[f]])
}
rownames(cd) <- cd$sample

# optional gene annotation + flags
gi <- if (!is.null(opt$gene_info) && file.exists(opt$gene_info)) fread(opt$gene_info) else data.table(gene_id=gene_ids)
flags <- if (!is.null(opt$flag_genes) && file.exists(opt$flag_genes)) fread(opt$flag_genes) else data.table(gene_id=character(), category=character())

cat(sprintf(">>> %d genes x %d samples; design %s\n", nrow(counts), ncol(counts), opt$design))

# ---- fit -------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(counts, cd, design = design_fm)
dds <- dds[rowSums(counts(dds)) >= opt$min_count, ]
dds <- DESeq(dds)
cat(">>> resultsNames:\n"); print(resultsNames(dds))
# vst() fits on >=1000 genes by default; fall back to the exact transform on
# small datasets (e.g. a downsampled test) where vst() errors.
vsd <- tryCatch(vst(dds, blind=TRUE),
                error=function(e) varianceStabilizingTransformation(dds, blind=TRUE))
saveRDS(dds, file.path(opt$outdir,"dds.rds")); saveRDS(vsd, file.path(opt$outdir,"vsd.rds"))

# ---- QC figures ------------------------------------------------------------
intg <- intersect(fct, names(cd)); intg <- head(intg, 3)
pc <- plotPCA(vsd, intgroup=intg, returnData=TRUE); pv <- round(100*attr(pc,"percentVar"))
aes_pca <- if (length(intg)>=2) aes(PC1,PC2,color=.data[[intg[1]]],shape=.data[[intg[2]]]) else aes(PC1,PC2,color=.data[[intg[1]]])
ggplot(pc, aes_pca) + geom_point(size=4) +
  labs(x=paste0("PC1 (",pv[1],"%)"), y=paste0("PC2 (",pv[2],"%)"), title="PCA (variance-stabilised counts)") + theme_bw()
ggsave(file.path(figdir,"PCA.png"), width=7, height=5, dpi=150)
sd <- dist(t(assay(vsd)))
pheatmap(as.matrix(sd), clustering_distance_rows=sd, clustering_distance_cols=sd,
         col=colorRampPalette(rev(brewer.pal(9,"Blues")))(255), main="Sample distance (VST)",
         filename=file.path(figdir,"sample_distance_heatmap.png"), width=7.3, height=6.3)
png(file.path(figdir,"dispersion.png"), width=900, height=750, res=150); plotDispEsts(dds); dev.off()

# ---- contrast engine -------------------------------------------------------
annotate <- function(res, shr, id) {
  res$log2FoldChange <- shr$log2FoldChange; res$lfcSE <- shr$lfcSE
  df <- as.data.frame(res); df$gene_id <- rownames(df)
  df <- merge(df, gi, by="gene_id", all.x=TRUE)
  for (cat in unique(flags$category)) df[[paste0("is_",cat)]] <- df$gene_id %in% flags[category==cat]$gene_id
  df <- df[order(df$padj), ]; fwrite(df, file.path(opt$outdir, paste0(id,".csv")))
  ns <- sum(df$padj < opt$padj, na.rm=TRUE)
  cat(sprintf(">>> %-28s padj<%.2g: %4d  (up %d / down %d)\n", id, opt$padj, ns,
      sum(df$padj<opt$padj & df$log2FoldChange>0, na.rm=TRUE), sum(df$padj<opt$padj & df$log2FoldChange<0, na.rm=TRUE)))
  png(file.path(figdir, paste0("MA_",id,".png")), width=900, height=750, res=150)
  DESeq2::plotMA(shr, ylim=c(-5,5), main=paste0("MA: ",id)); dev.off()
  v <- df[!is.na(df$padj),]; v$sig <- v$padj < opt$padj & abs(v$log2FoldChange) > 1
  lab <- head(v[order(v$padj),], 15); nm <- if ("std_name" %in% names(v)) "std_name" else "gene_id"
  ggplot(v, aes(log2FoldChange, -log10(padj), color=sig)) + geom_point(alpha=.5, size=1.2) +
    scale_color_manual(values=c(`FALSE`="grey70",`TRUE`="firebrick")) +
    geom_text_repel(data=lab, aes(label=.data[[nm]]), size=3, max.overlaps=20, show.legend=FALSE) +
    geom_hline(yintercept=-log10(opt$padj), linetype=2) + geom_vline(xintercept=c(-1,1), linetype=2) +
    labs(title=id, x="shrunken log2FC", y="-log10 padj") + theme_bw()
  ggsave(file.path(figdir, paste0("volcano_",id,".png")), width=8, height=6, dpi=150)
  list(df=df, nsig=ns)
}
ct <- fread(opt$contrasts); res_list <- list(); summ <- list()
for (i in seq_len(nrow(ct))) {
  r <- ct[i]; id <- r$id
  if (r$type == "pairwise") {
    co <- c(r$factor, r$target, r$reference)
    res <- results(dds, contrast=co); shr <- lfcShrink(dds, contrast=co, type="ashr", quiet=TRUE)
  } else if (r$type == "simple") {
    wf <- strsplit(r$within, "=")[[1]]
    d2 <- dds; d2[[wf[1]]] <- relevel(d2[[wf[1]]], ref=wf[2]); d2[[r$factor]] <- relevel(d2[[r$factor]], ref=r$reference)
    d2 <- DESeq(d2, quiet=TRUE); coef <- paste0(r$factor,"_",r$target,"_vs_",r$reference)
    if (!coef %in% resultsNames(d2)) stop("coef not found for ",id,": ",coef)
    res <- results(d2, name=coef); shr <- lfcShrink(d2, coef=coef, type="apeglm", quiet=TRUE)
  } else if (r$type == "interaction") {
    fs <- strsplit(r$factor, ":")[[1]]; ts <- strsplit(r$target, ":")[[1]]
    coef <- paste0(fs[1], ts[1], ".", fs[2], ts[2])
    if (!coef %in% resultsNames(dds)) stop("interaction coef not found for ",id,": ",coef,
        " (have: ", paste(resultsNames(dds),collapse=", "),")")
    res <- results(dds, name=coef); shr <- lfcShrink(dds, coef=coef, type="apeglm", quiet=TRUE)
  } else stop("unknown contrast type: ", r$type)
  a <- annotate(res, shr, id); res_list[[id]] <- a$df
  summ[[id]] <- data.table(contrast=id, type=r$type, n_sig=a$nsig)
}
saveRDS(res_list, file.path(opt$outdir,"res_list.rds"))
fwrite(rbindlist(summ), file.path(opt$outdir,"deg_summary.csv"))
cat(">>> DESeq2 complete.\n")
