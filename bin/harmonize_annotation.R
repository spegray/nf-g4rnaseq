#!/usr/bin/env Rscript
# =============================================================================
# harmonize_annotation.R  -  Reconcile genome + annotation (+ BAM) seqnames
# Generalises analysis/scripts/01_prepare_annotation.sh. Makes the genome FASTA
# and GFF use the SAME chromosome names as the reads' coordinate system, then
# emits the featureCounts SAF, a gene_info table, and an indexed genome.
#   --harmonize length : auto-match contigs by LENGTH to the target naming
#                        (target = BAM @SQ names if --bam given, else genome names)
#   --harmonize map    : use an explicit 2-column old<TAB>new map (--map)
#   --harmonize none   : assume names already match (just build SAF/gene_info/fa)
# =============================================================================
suppressPackageStartupMessages({ library(optparse); library(Biostrings); library(data.table) })
opt <- parse_args(OptionParser(option_list=list(
  make_option("--gff",   type="character", help="GFF3 annotation (may be gzipped; may embed ##FASTA)"),
  make_option("--fasta", type="character", help="genome FASTA (may be gzipped)"),
  make_option("--bam",   type="character", default=NULL, help="optional BAM whose @SQ names are the target"),
  make_option("--samtools", type="character", default="samtools", help="samtools executable"),
  make_option("--harmonize", type="character", default="length", help="length | map | none"),
  make_option("--map",   type="character", default=NULL, help="2-col TSV old<TAB>new (for --harmonize map)"),
  make_option("--outdir",type="character", default=".")
)))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

genome <- readDNAStringSet(opt$fasta); names(genome) <- sub("\\s.*","",names(genome))
glen <- setNames(as.integer(width(genome)), names(genome))

# target naming (name -> length)
if (!is.null(opt$bam)) {
  hdr <- system(sprintf("%s view -H %s", opt$samtools, shQuote(opt$bam)), intern=TRUE)
  sq <- hdr[grepl("^@SQ", hdr)]
  tname <- sub(".*SN:([^\t]+).*", "\\1", sq); tlen <- as.integer(sub(".*LN:([0-9]+).*", "\\1", sq))
  target <- setNames(tlen, tname)
} else target <- glen

# build a function: given a vector of (name->length), map names to target by chosen mode
make_map <- function(src) {
  if (opt$harmonize == "none") return(setNames(names(src), names(src)))
  if (opt$harmonize == "map") { m <- fread(opt$map, header=FALSE); return(setNames(m$V2, m$V1)) }
  # length mode: match each src contig to the unique target contig of equal length
  tl <- target; out <- setNames(names(src), names(src))
  for (nm in names(src)) { hit <- names(tl)[tl == src[[nm]]]; if (length(hit)==1) out[[nm]] <- hit }
  out
}
gmap <- make_map(glen)
names(genome) <- unname(gmap[names(genome)])
writeXStringSet(genome, file.path(opt$outdir,"genome.fa"))
system(sprintf("%s faidx %s", opt$samtools, shQuote(file.path(opt$outdir,"genome.fa"))))

# ---- read GFF annotation rows (stop at ##FASTA), get per-seqid length --------
con <- if (grepl("\\.gz$", opt$gff)) gzfile(opt$gff) else file(opt$gff)
gff <- readLines(con); close(con)
fasta_start <- which(grepl("^##FASTA", gff)); if (length(fasta_start)) gff <- gff[seq_len(fasta_start[1]-1)]
body <- gff[!grepl("^#", gff) & nchar(gff)>0]
dt <- fread(text=paste(body, collapse="\n"), header=FALSE, sep="\t", quote="")
setnames(dt, paste0("V",1:9)[1:ncol(dt)])
slen <- dt[, .(len=max(V5)), by=V1]; src <- setNames(slen$len, slen$V1)
amap <- make_map(src)
dt[, V1 := amap[V1]]
gff_out <- file.path(opt$outdir,"annotation.gff3")
fwrite(dt, gff_out, sep="\t", col.names=FALSE, quote=FALSE)

# ---- SAF + gene_info from gene-class features --------------------------------
g <- dt[grepl("gene$", V3)]
getattr <- function(a, key) { m <- regmatches(a, regexpr(sprintf("%s=[^;]+", key), a)); ifelse(length(m)>0 && nchar(m), sub(".*=","",m), NA_character_) }
g[, gene_id := vapply(V9, getattr, "", key="ID")]
g[, std_name := vapply(V9, getattr, "", key="gene")]
g[, orf_class := vapply(V9, getattr, "", key="orf_classification")]
g <- g[!is.na(gene_id)]; g[is.na(std_name), std_name := gene_id]; g[is.na(orf_class), orf_class := "NA"]
fwrite(g[, .(GeneID=gene_id, Chr=V1, Start=V4, End=V5, Strand=V7)], file.path(opt$outdir,"genes.saf"), sep="\t")
fwrite(g[, .(gene_id, std_name, orf_class, chr=V1, start=V4, end=V5, strand=V7, biotype=V3)],
       file.path(opt$outdir,"gene_info.tsv"), sep="\t")
cat(sprintf(">>> harmonized to %d contigs; %d gene loci; SAF + gene_info + genome.fa written to %s\n",
    length(genome), nrow(g), opt$outdir))
