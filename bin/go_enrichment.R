#!/usr/bin/env Rscript
# =============================================================================
# go_enrichment.R  -  GO + KEGG over-representation per contrast (up & down)
# Generalises analysis/scripts/09 + 10. Organism via --orgdb / --kegg_organism.
# Skips gracefully if --orgdb is not installed (so the rest of the pipeline runs).
# =============================================================================
suppressPackageStartupMessages({ library(optparse); library(data.table) })
opt <- parse_args(OptionParser(option_list=list(
  make_option("--deseq2_dir", type="character"),
  make_option("--orgdb",      type="character", default=NULL, help="e.g. org.Sc.sgd.db, org.Hs.eg.db"),
  make_option("--kegg_organism", type="character", default=NULL, help="e.g. sce, hsa (KEGG code; needs internet)"),
  make_option("--keytype",    type="character", default="ORF"),
  make_option("--flag_genes", type="character", default=NULL, help="genes to EXCLUDE (confounders)"),
  make_option("--ontologies", type="character", default="BP,MF,CC"),
  make_option("--padj",       type="double", default=0.05),
  make_option("--outdir",     type="character", default=".")
)))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE); figdir<-file.path(opt$outdir,"figures"); dir.create(figdir,showWarnings=FALSE)
if (is.null(opt$orgdb) || !requireNamespace(opt$orgdb, quietly=TRUE)) {
  cat(">>> No usable --orgdb ('", opt$orgdb, "') -> skipping GO/KEGG (rest of pipeline unaffected).\n", sep=""); quit(save="no", status=0)
}
suppressPackageStartupMessages({ library(clusterProfiler); library(opt$orgdb, character.only=TRUE); library(ggplot2) })
orgdb <- get(opt$orgdb)
rl <- readRDS(file.path(opt$deseq2_dir,"res_list.rds"))
flag_ids <- if (!is.null(opt$flag_genes) && file.exists(opt$flag_genes)) fread(opt$flag_genes)$gene_id else character()
onts <- strsplit(opt$ontologies, ",")[[1]]

run_go <- function(genes, universe, tag, ont) {
  if (length(genes)<5) return(invisible())
  ego <- tryCatch(enrichGO(genes, orgdb, keyType=opt$keytype, ont=ont, universe=universe,
            pAdjustMethod="BH", pvalueCutoff=opt$padj, qvalueCutoff=0.2), error=function(e) NULL)
  if (is.null(ego) || nrow(as.data.frame(ego))==0) { cat(sprintf("   GO-%s %-28s: none\n",ont,tag)); return(invisible()) }
  if (ont=="BP") ego <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
  df <- as.data.frame(ego); fwrite(df, file.path(opt$outdir, sprintf("GO_%s_%s.csv",ont,tag)))
  ggsave(file.path(figdir, sprintf("dotplot_GO%s_%s.png",ont,tag)),
         clusterProfiler::dotplot(ego, showCategory=min(15,nrow(df)))+ggtitle(sprintf("GO:%s %s",ont,tag)), width=9,height=6,dpi=150)
  cat(sprintf("   GO-%s %-28s: %d terms | top: %s (padj=%.1e)\n",ont,tag,nrow(df),df$Description[1],df$p.adjust[1]))
}
for (cn in names(rl)) {
  d <- as.data.table(rl[[cn]])[!is.na(padj)]; uni <- d$gene_id
  sig <- d[padj<opt$padj & !gene_id %in% flag_ids]
  for (ont in onts) {
    run_go(sig[log2FoldChange<0]$gene_id, uni, paste0(cn,"_down"), ont)
    run_go(sig[log2FoldChange>0]$gene_id, uni, paste0(cn,"_up"),   ont)
  }
  if (!is.null(opt$kegg_organism)) {
    kk <- tryCatch(enrichKEGG(sig$gene_id, organism=opt$kegg_organism, universe=uni, pvalueCutoff=0.1), error=function(e) NULL)
    if (!is.null(kk) && nrow(as.data.frame(kk))>0) { fwrite(as.data.frame(kk), file.path(opt$outdir, sprintf("KEGG_%s.csv",cn)))
      cat(sprintf("   KEGG %-28s: %d pathways | top: %s\n", cn, nrow(as.data.frame(kk)), as.data.frame(kk)$Description[1])) }
  }
}
cat(">>> go_enrichment complete.\n")
