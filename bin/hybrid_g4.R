#!/usr/bin/env Rscript
# =============================================================================
# hybrid_g4.R  -  Co-transcriptional RNA:DNA hybrid-G4 (PHQS) scan
# Generalises analysis/scripts/07_hybrid_g4.R. Implements the PHQS model of
# Zheng et al. 2013 (NAR): count clustered G-tracts on the SENSE (coding) strand.
#   >= min_tracts G-runs (each >= run_len G) within loop_max gaps, span <= span_max
#   -> hybrid-capable. 2-3 tracts (cannot fold a 4-tract DNA-only G4) -> hybrid-specific.
# Organism-agnostic. Optional expression table flags which loci are transcribed.
# =============================================================================
suppressPackageStartupMessages({ library(optparse); library(Biostrings); library(data.table) })
opt <- parse_args(OptionParser(option_list=list(
  make_option("--fasta",     type="character", help="genome FASTA"),
  make_option("--gene_info", type="character", help="TSV: gene_id, std_name, biotype, chr, start, end, strand"),
  make_option("--expr",      type="character", default=NULL, help="optional TSV: gene_id, baseMean (to flag transcribed loci)"),
  make_option("--min_expr",  type="double",  default=10, help="baseMean threshold for 'expressed'"),
  make_option("--run_len",   type="integer", default=3L,  help="min consecutive G's per tract"),
  make_option("--loop_max",  type="integer", default=12L, help="max gap (nt) between consecutive tracts in a cluster"),
  make_option("--span_max",  type="integer", default=100L,help="max nt a foldable cluster may span"),
  make_option("--min_tracts",type="integer", default=2L,  help="min tracts to be hybrid-capable (2 = 2 DNA:2 RNA)"),
  make_option("--outdir",    type="character", default=".")
)))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
genome <- readDNAStringSet(opt$fasta); names(genome) <- sub(" .*","",names(genome))
gi <- fread(opt$gene_info); gi <- gi[chr %in% names(genome)]
if (!is.null(opt$expr) && file.exists(opt$expr)) {
  e <- fread(opt$expr); gi[, baseMean := e$baseMean[match(gene_id, e$gene_id)]]
} else gi[, baseMean := NA_real_]
runpat <- sprintf("G{%d,}", opt$run_len)

scan_gene <- function(seq) {
  m <- gregexpr(runpat, seq, perl=TRUE)[[1]]
  gfrac <- letterFrequency(DNAString(seq),"G",as.prob=TRUE)[[1]]
  if (m[1]==-1L) return(list(nt=0L,nhyb=0L,nspec=0L,ncan=0L,maxk=0L,gfrac=gfrac))
  s <- as.integer(m); e <- s + attr(m,"match.length") - 1L; cl <- integer(length(s)); cl[1] <- 1L
  if (length(s)>1) for (i in 2:length(s)) cl[i] <- if ((s[i]-e[i-1]-1L)>opt$loop_max) cl[i-1]+1L else cl[i-1]
  dt <- data.table(s,e,cl)[, .(k=.N, span=max(e)-min(s)+1L), by=cl]
  f <- dt[k>=opt$min_tracts & span<=opt$span_max]
  list(nt=length(s), nhyb=nrow(f), nspec=nrow(f[k%in%2:3]), ncan=nrow(f[k>=4]),
       maxk=if(nrow(f)) max(f$k) else 0L, gfrac=gfrac)
}
out <- vector("list", nrow(gi))
for (i in seq_len(nrow(gi))) {
  g <- gi[i]; sq <- subseq(genome[[g$chr]], g$start, g$end); if (g$strand=="-") sq <- reverseComplement(sq)
  r <- scan_gene(as.character(sq))
  out[[i]] <- data.table(gene_id=g$gene_id, std_name=g$std_name, biotype=g$biotype, strand=g$strand,
    length=g$end-g$start+1L, baseMean=g$baseMean, sense_G_frac=round(r$gfrac,3),
    sense_G3_tracts=r$nt, hybrid_capable_clusters=r$nhyb, hybrid_specific_clusters=r$nspec,
    canonical_sense_clusters=r$ncan, max_tracts_in_cluster=r$maxk)
}
H <- rbindlist(out)
H[, hybrid_capable  := hybrid_capable_clusters>0]
H[, hybrid_specific := hybrid_specific_clusters>0 & canonical_sense_clusters==0]
H[, expressed := !is.na(baseMean) & baseMean>=opt$min_expr]
fwrite(H, file.path(opt$outdir,"gene_hybrid_g4.tsv"), sep="\t")
fwrite(H[hybrid_specific==TRUE][order(sense_G_frac)], file.path(opt$outdir,"hybrid_specific_genes.tsv"), sep="\t")
cat(sprintf(">>> hybrid-capable: %d (%.1f%%); hybrid-specific: %d (%.1f%%); expressed hybrid-specific: %d\n",
    sum(H$hybrid_capable), 100*mean(H$hybrid_capable), sum(H$hybrid_specific), 100*mean(H$hybrid_specific),
    sum(H$hybrid_specific & H$expressed, na.rm=TRUE)))
cat(">>> hybrid_g4 complete.\n")
