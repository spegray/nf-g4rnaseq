// FASTQ entry path: QC/trim -> STAR align -> index. (UMI handling via umi_tools
// can be inserted between FASTP and STAR for UMI libraries; see docs/usage.md.)
// These use standard biocontainers and run under -profile docker/singularity.

process FASTP {
    label 'process_medium'
    container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_2'
    input:  tuple val(sample), path(reads)
    output: tuple val(sample), path("${sample}.trim*.fastq.gz"), emit: reads
            path "${sample}.fastp.json", emit: json
    script:
      if (reads instanceof List && reads.size() == 2)
        """
        fastp -i ${reads[0]} -I ${reads[1]} -o ${sample}.trim_1.fastq.gz -O ${sample}.trim_2.fastq.gz \\
              -j ${sample}.fastp.json -w ${task.cpus}
        """
      else
        """
        fastp -i ${reads} -o ${sample}.trim.fastq.gz -j ${sample}.fastp.json -w ${task.cpus}
        """
}

process STAR_INDEX {
    label 'process_high'
    container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_1'
    input:  path fasta
    output: path "star_index", emit: index
    script:
      """
      mkdir star_index
      STAR --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles ${fasta} \\
           --runThreadN ${task.cpus} --genomeSAindexNbases 11
      """
}

process STAR_ALIGN {
    label 'process_high'
    container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_1'
    input:  tuple val(sample), path(reads)
            path index
    output: tuple val(sample), path("${sample}.Aligned.sortedByCoord.out.bam"), emit: bam
            path "${sample}.Log.final.out", emit: log
    script:
      """
      STAR --genomeDir ${index} --readFilesIn ${reads} --readFilesCommand zcat --runThreadN ${task.cpus} \\
           --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${sample}. \\
           --outFilterIntronMotifs RemoveNoncanonical --outSAMattributes NH HI AS nM
      """
}

process SAMTOOLS_INDEX {
    label 'process_low'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'
    input:  tuple val(sample), path(bam)
    output: tuple val(sample), path(bam), emit: bam
    script: "samtools index ${bam}"
}
