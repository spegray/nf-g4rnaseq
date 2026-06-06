#!/usr/bin/env Rscript
# =============================================================================
# g4_enrichment.R  -  Test whether DEGs are enriched for proximal/hybrid G4s
# Generalises analysis/scripts/06 + 11 + 12 (directional + ribosome-excluded).
# For every contrast: Fisher test of has_g4_highconf / hybrid_capable /
# hybrid_specific among DEGs (all/up/down) vs the tested-gene background, with
# flagged (e.g. mating-type/construction) genes excluded. If --orgdb is given,
# also reports the result with a GO-defined gene set (default: ribosome/
# translation) removed, as a robustness check.
# =============================================================================
suppressPackageStartupMessages({ library(optparse); library(data.table) })
opt <- parse_args(OptionParser(option_list=list(
  make_option("--deseq2_dir",  type="character", help="dir with res_list.rds"),
  make_option("--g4_features", type="character"),
  make_option("--hybrid",      type="character"),
  make_option("--flag_genes",  type="character", default=NULL),
  make_option("--orgdb",       type="character", default=NULL, help="OrgDb for the robustness gene set (optional)"),
  make_option("--robustness_go", type="character", default="GO:0022613,GO:0042254,GO:0006364,GO:0006412,GO:0003735",
              help="GO IDs whose genes are removed in the robustness check"),
  make_option("--keytype",     type="character", default="ORF"),
  make_option("--padj",        type="double", default=0.05),
  make_option("--outdir",      type="character", default=".")
)))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
rl <- readRDS(file.path(opt$deseq2_dir,"res_list.rds"))
g4 <- fread(opt$g4_features); hy <- fread(opt$hybrid)
feat <- merge(g4[, .(gene_id, has_g4_highconf, n_g4_coding)],
              hy[, .(gene_id, hybrid_capable, hybrid_specific)], by="gene_id", all=TRUE)
feat[, coding_g4 := n_g4_coding>0]
flag_ids <- if (!is.null(opt$flag_genes) && file.exists(opt$flag_genes)) fread(opt$flag_genes)$gene_id else character()

# optional robustness gene set
ribo <- character()
if (!is.null(opt$orgdb) && requireNamespace(opt$orgdb, quietly=TRUE)) {
  suppressPackageStartupMessages(library(opt$orgdb, character.only=TRUE)); library(AnnotationDbi)
  ribo <- tryCatch(unique(AnnotationDbi::select(get(opt$orgdb), keys=strsplit(opt$robustness_go,",")[[1]],
            keytype="GOALL", columns=opt$keytype)[[opt$keytype]]), error=function(e) character())
  ribo <- ribo[!is.na(ribo)]
}
flags <- c("has_g4_highconf","hybrid_capable","hybrid_specific","coding_g4")
fish <- function(deg, bg, f) {
  a<-sum(deg[[f]],na.rm=TRUE); b<-nrow(deg)-a; c<-sum(bg[[f]],na.rm=TRUE); dd<-nrow(bg)-c
  ft<-fisher.test(matrix(c(a,b,c,dd),2))
  data.table(flag=f, pct_deg=round(100*a/nrow(deg),1), pct_bg=round(100*c/nrow(bg),1),
             odds_ratio=round(unname(ft$estimate),2), p=signif(ft$p.value,3), n_deg=nrow(deg))
}
res <- list()
for (cn in names(rl)) {
  d <- as.data.table(rl[[cn]]); d <- merge(d[!is.na(padj)], feat, by="gene_id", all.x=TRUE)
  flagcol <- grep("^is_", names(d), value=TRUE)
  keep <- if (length(flagcol)) rowSums(as.matrix(d[, ..flagcol]))==0 else rep(TRUE,nrow(d))
  uni <- d[keep]; sig <- uni[padj<opt$padj]
  for (subset_name in c("with_ribo","no_ribo")) {
    u <- if (subset_name=="no_ribo") uni[!gene_id %in% ribo] else uni
    s <- if (subset_name=="no_ribo") sig[!gene_id %in% ribo] else sig
    if (subset_name=="no_ribo" && length(ribo)==0) next
    for (dir in c("all","up","down")) {
      gg <- if (dir=="all") s else if (dir=="up") s[log2FoldChange>0] else s[log2FoldChange<0]
      if (nrow(gg)<10) next
      bg <- u[!gene_id %in% gg$gene_id]
      for (f in flags) res[[length(res)+1]] <- cbind(contrast=cn, set=subset_name, direction=dir, fish(gg,bg,f))
    }
  }
}
enr <- rbindlist(res, fill=TRUE)
fwrite(enr, file.path(opt$outdir,"g4_deg_enrichment.csv"))
cat(">>> g4_deg_enrichment.csv:", nrow(enr), "rows. Key (high-conf G4, all DEGs, with ribo):\n")
print(enr[set=="with_ribo" & direction=="all" & flag=="has_g4_highconf"][order(contrast)])
cat(">>> g4_enrichment complete.\n")
