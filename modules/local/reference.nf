// Reference preparation: harmonize genome + annotation (+ optional BAM) seqnames,
// then emit indexed genome, gene_info, featureCounts SAF, and renamed GFF.
process HARMONIZE {
    label 'process_low'
    container params.r_container          // needs Biostrings + samtools on PATH
    input:
      path fasta
      path gff
      path bam            // assets/NO_FILE if none (counts/no-target entry)
      path seqmap         // assets/NO_FILE unless harmonize=map
    output:
      path "reference/genome.fa",        emit: fasta
      path "reference/genome.fa.fai",    emit: fai
      path "reference/gene_info.tsv",    emit: gene_info
      path "reference/genes.saf",        emit: saf
      path "reference/annotation.gff3",  emit: gff
    script:
      def bamopt = bam.name    != 'NO_BAM' ? "--bam ${bam}"   : ''
      def mapopt = seqmap.name != 'NO_MAP' ? "--map ${seqmap}": ''
      """
      mkdir -p reference
      harmonize_annotation.R --gff ${gff} --fasta ${fasta} ${bamopt} \\
        --harmonize ${params.harmonize} ${mapopt} --samtools samtools --outdir reference
      """
}
