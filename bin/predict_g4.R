#!/usr/bin/env Rscript
# =============================================================================
# predict_g4.R  -  Genome-wide G4 prediction + per-gene proximal-G4 features
# Generalises analysis/scripts/05_predict_g4.R (+ the feature part of 06).
# Organism-agnostic: needs only a genome FASTA and a gene_info table.
# =============================================================================
suppressPackageStartupMessages({
  library(optparse); library(pqsfinder); library(Biostrings); library(GenomicRanges); library(data.table)
})
opt <- parse_args(OptionParser(option_list=list(
  make_option("--fasta",     type="character", help="genome FASTA (seqnames matching the annotation/BAMs)"),
  make_option("--gene_info", type="character", default=NULL, help="TSV: gene_id, chr, start, end, strand (for per-gene features)"),
  make_option("--min_score", type="integer",   default=25L, help="pqsfinder min score (permissive=25, default=47)"),
  make_option("--promoter",  type="integer",   default=500L, help="bp upstream of TSS included as 'proximal'"),
  make_option("--highconf_score", type="integer", default=47L, help="score threshold for the high-confidence subset"),
  make_option("--outdir",    type="character", default=".")
)))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
genome <- readDNAStringSet(opt$fasta); names(genome) <- sub(" .*","",names(genome))
chrlen <- setNames(width(genome), names(genome))

# ---- pqsfinder on both strands, per chromosome -----------------------------
pqs <- list()
for (chr in names(genome)) {
  pv <- pqsfinder(genome[[chr]], strand="*", min_score=opt$min_score); g <- as(pv,"GRanges")
  g2 <- GRanges(chr, ranges(g), strand(g)); mcols(g2) <- mcols(g); seqlevels(g2) <- names(genome)
  pqs[[chr]] <- g2
}
g4 <- sort(do.call(c, unname(pqs))); seqlevels(g4) <- names(genome); seqlengths(g4) <- chrlen
mcols(g4)$g4_id <- sprintf("G4_%06d", seq_along(g4))
saveRDS(g4, file.path(opt$outdir,"g4_sites.rds"))
fwrite(as.data.table(g4)[,.(g4_id,seqnames,start,end,strand,score,n_tetrads=nt,bulges=nb,mismatches=nm)],
       file.path(opt$outdir,"g4_sites.tsv"), sep="\t")
fwrite(as.data.table(g4)[,.(seqnames,start=start-1L,end,name=g4_id,score=pmin(1000L,as.integer(score)),strand)],
       file.path(opt$outdir,"g4_sites.bed"), sep="\t", col.names=FALSE)
cat(sprintf(">>> pqsfinder G4s (score>=%d): %d (+:%d -:%d); canonical(>=3 tetrad):%d; two-quartet:%d; high-conf(>=%d):%d\n",
    opt$min_score, length(g4), sum(strand(g4)=="+"), sum(strand(g4)=="-"),
    sum(mcols(g4)$nt>=3), sum(mcols(g4)$nt==2), opt$highconf_score, sum(mcols(g4)$score>=opt$highconf_score)))

# ---- per-gene proximal-G4 features (promoter + gene body, strand-aware) -----
if (!is.null(opt$gene_info) && file.exists(opt$gene_info)) {
  gi <- fread(opt$gene_info); gi <- gi[chr %in% names(chrlen)]
  genes <- GRanges(gi$chr, IRanges(gi$start, gi$end), strand=gi$strand, gene_id=gi$gene_id)
  seqlevels(genes) <- names(chrlen); seqlengths(genes) <- chrlen
  st <- as.character(strand(genes)); gs <- start(genes); ge <- end(genes)
  rs <- ifelse(st=="+", pmax(1L, gs-opt$promoter), gs)
  re <- ifelse(st=="+", ge, pmin(chrlen[as.character(seqnames(genes))], ge+opt$promoter))
  region <- GRanges(seqnames(genes), IRanges(rs,re), strand=strand(genes))
  count_g4 <- function(sub) {
    ov <- findOverlaps(region, sub, ignore.strand=TRUE)
    dt <- data.table(gene=queryHits(ov), gstr=st[queryHits(ov)], g4str=as.character(strand(sub))[subjectHits(ov)])
    dt[, rel := ifelse(g4str==gstr,"coding","template")]
    list(n=tabulate(dt$gene,nbins=length(genes)),
         cod=tabulate(dt[rel=="coding"]$gene,nbins=length(genes)),
         tem=tabulate(dt[rel=="template"]$gene,nbins=length(genes)))
  }
  a <- count_g4(g4); cano <- count_g4(g4[mcols(g4)$nt>=3]); hi <- count_g4(g4[mcols(g4)$score>=opt$highconf_score])
  feat <- data.table(gene_id=genes$gene_id, n_g4=a$n, n_g4_coding=a$cod, n_g4_template=a$tem,
                     n_g4_canonical=cano$n, n_g4_highconf=hi$n,
                     has_g4=a$n>0, has_g4_canonical=cano$n>0, has_g4_highconf=hi$n>0)
  fwrite(feat, file.path(opt$outdir,"gene_g4_features.tsv"), sep="\t")
  cat(sprintf(">>> per-gene G4 features: %.1f%% genes have a proximal G4 (any); %.1f%% high-conf\n",
      100*mean(feat$has_g4), 100*mean(feat$has_g4_highconf)))
}
cat(">>> predict_g4 complete.\n")
