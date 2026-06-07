// featureCounts quantification from BAMs. Strandedness is auto-detected
// (testing -s 0/1/2 on the first BAM and choosing the highest assignment) unless
// params.strandedness is set to 0/1/2. Count columns are renamed to the sample
// IDs (BAMs are passed in sample-sheet order).
process FEATURECOUNTS {
    label 'process_high'
    container 'quay.io/biocontainers/subread:2.0.6--he4a0461_2'
    input:
      path bams           // ordered to match `samples`
      val  samples        // ordered list of sample IDs
      path saf
    output:
      path "counts.tsv",                emit: counts
      path "gene_counts.txt.summary",   emit: summary
    script:
      def names = samples.join(',')
      """
      S=${params.strandedness}
      if [ "\$S" = "auto" ]; then
        best=0; bestn=-1
        for s in 0 1 2; do
          featureCounts -F SAF -a ${saf} -s \$s -T ${task.cpus} -o probe_\$s.txt \$(echo ${bams} | tr ' ' '\\n' | head -1) >/dev/null 2>&1 || true
          n=\$(awk '\$1=="Assigned"{print \$2}' probe_\$s.txt.summary 2>/dev/null || echo 0)
          if [ "\$n" -gt "\$bestn" ]; then bestn=\$n; best=\$s; fi
        done; S=\$best
        echo ">>> auto-detected strandedness: -s \$S"
      fi
      featureCounts -F SAF -a ${saf} -s \$S -T ${task.cpus} -o gene_counts.txt ${bams}
      Rscript -e 'library(data.table); fc<-fread("gene_counts.txt",skip=1); m<-fc[,c(1,7:ncol(fc)),with=FALSE]; nm<-strsplit("${names}",",")[[1]]; stopifnot((ncol(fc)-6)==length(nm)); setnames(m, c("gene_id", nm)); fwrite(m,"counts.tsv",sep="\\t")'
      """
}
