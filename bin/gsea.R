#!/usr/bin/env Rscript
# =============================================================================
# gsea.R  -  Threshold-free GSEA of every contrast (genes ranked by Wald stat)
# Generalises analysis/scripts/13 + 13b + 14. Runs uniformly on within-factor,
# pairwise AND interaction contrasts (they are all rows of res_list). For each:
#   (1) gseGO BP (if --orgdb available)
#   (2) custom GSEA of G4 gene sets (hybrid_specific/capable, high-conf proximal,
#       coding-strand) -> does the contrast coordinately shift G4-bearing genes?
# Plus a summary NES bar chart of the G4 sets across contrasts.
# =============================================================================
suppressPackageStartupMessages({ library(optparse); library(clusterProfiler); library(limma); library(ggplot2); library(data.table) })
opt <- parse_args(OptionParser(option_list=list(
  make_option("--deseq2_dir",  type="character"),
  make_option("--g4_features", type="character"),
  make_option("--hybrid",      type="character"),
  make_option("--orgdb",       type="character", default=NULL),
  make_option("--keytype",     type="character", default="ORF"),
  make_option("--outdir",      type="character", default=".")
)))
set.seed(1); dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE); figdir<-file.path(opt$outdir,"figures"); dir.create(figdir,showWarnings=FALSE)
rl <- readRDS(file.path(opt$deseq2_dir,"res_list.rds"))
g4 <- fread(opt$g4_features); hy <- fread(opt$hybrid)
t2g <- rbindlist(list(
  data.table(term="hybrid_specific",      gene=hy[hybrid_specific==TRUE]$gene_id),
  data.table(term="hybrid_capable",       gene=hy[hybrid_capable==TRUE]$gene_id),
  data.table(term="proximal_highconf_G4", gene=g4[has_g4_highconf==TRUE]$gene_id),
  data.table(term="coding_strand_G4",     gene=g4[n_g4_coding>0]$gene_id)))
have_org <- !is.null(opt$orgdb) && requireNamespace(opt$orgdb, quietly=TRUE)
if (have_org) suppressPackageStartupMessages(library(opt$orgdb, character.only=TRUE))
g4collect <- list()
for (cn in names(rl)) {
  d <- as.data.table(rl[[cn]])[!is.na(stat)]; gl <- d$stat; names(gl) <- d$gene_id; gl <- sort(gl, decreasing=TRUE)
  if (have_org) {
    eg <- tryCatch(gseGO(gl, OrgDb=get(opt$orgdb), keyType=opt$keytype, ont="BP", minGSSize=10, maxGSSize=500,
                         pvalueCutoff=0.1, eps=0, verbose=FALSE), error=function(e) NULL)
    if (!is.null(eg) && nrow(as.data.frame(eg))>0) fwrite(as.data.frame(eg)[,c("ID","Description","setSize","NES","pvalue","p.adjust")],
                                                          file.path(opt$outdir, sprintf("gseGO_%s.csv",cn)))
  }
  # G4 gene sets are large (often thousands of genes), where fgsea's NES/p are
  # inflated by set size and ignore inter-gene correlation. Use limma::cameraPR
  # instead: a COMPETITIVE, correlation-aware test on the ranked statistic,
  # calibrated for any set size (inter.gene.cor=0.01 is the recommended robust
  # default; use.ranks=TRUE = non-parametric). Reports Direction + BH-FDR; we add
  # the median in-set statistic as an interpretable effect size.
  idx <- lapply(split(t2g$gene, t2g$term), function(s) which(names(gl) %in% s))
  idx <- idx[sapply(idx, length) >= 10]
  if (length(idx)) {
    cam <- limma::cameraPR(gl, idx, use.ranks=TRUE, inter.gene.cor=0.01)
    cam$set <- rownames(cam)
    cam$median_stat <- sapply(cam$set, function(s) round(median(gl[idx[[s]]]), 3))
    out <- as.data.table(cam)[, .(ID=set, NGenes, Direction, median_stat,
                                  p=signif(PValue,3), padj=signif(FDR,3))]
    fwrite(out, file.path(opt$outdir, sprintf("g4_geneset_camera_%s.csv",cn)))
    g4collect[[cn]] <- cbind(contrast=cn, out[, .(ID, median_stat, Direction, padj)])
    cat(sprintf(">>> %-26s G4 cameraPR: %s\n", cn,
        paste(sprintf("%s %s(med=%.2f,padj=%.2g)", out$ID, out$Direction, out$median_stat, out$padj), collapse="  ")))
  }
}
if (length(g4collect)) {
  all <- rbindlist(g4collect); all[, sig := ifelse(padj<0.001,"***",ifelse(padj<0.01,"**",ifelse(padj<0.05,"*","")))]
  ggplot(all, aes(ID, median_stat, fill=median_stat)) + geom_col() +
    geom_text(aes(label=sig, vjust=ifelse(median_stat<0,1.2,-0.4))) +
    scale_fill_gradient2(low="firebrick", mid="grey90", high="steelblue", midpoint=0) +
    facet_wrap(~contrast) + coord_flip() +
    labs(title="G4 gene sets: competitive test (limma::cameraPR, correlation-aware)",
         subtitle="bar = median Wald stat in set (<0 = set down-regulated);  */**/*** = BH-FDR<0.05/0.01/0.001",
         x=NULL, y="median Wald statistic in set") +
    theme_bw() + theme(legend.position="none")
  ggsave(file.path(figdir,"G4_geneset_camera_byContrast.png"), width=11, height=7, dpi=150)
}
cat(">>> gsea complete.\n")
