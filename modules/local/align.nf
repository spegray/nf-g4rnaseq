// FASTQ entry path: QC/trim -> STAR align -> index. (UMI handling via umi_tools
// can be inserted between FASTP and STAR for UMI libraries; see docs/usage.md.)
// These use standard biocontainers and run under -profile docker/singularity.
// NOTE: STAR is a Linux-first tool; run the FASTQ path on Linux / in the
// container (the conda macOS-arm64 STAR build can silently read 0 input reads).

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
      # Decompress first and feed plain FASTQ. STAR's --readFilesCommand fails to
      # spawn on some platforms (e.g. macOS); decompressing is portable everywhere.
      R=""
      for f in ${reads}; do
        case "\$f" in
          *.gz) b=\$(basename "\$f" .gz); gunzip -c "\$f" > "\$b"; R="\$R \$b" ;;
          *)    R="\$R \$f" ;;
        esac
      done
      STAR --genomeDir ${index} --readFilesIn \$R --runThreadN ${task.cpus} \\
           --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${sample}. \\
           --outFilterIntronMotifs RemoveNoncanonical --outSAMattributes NH HI AS nM
      """
}

process SAMTOOLS_INDEX {
    label 'process_low'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'
    input:  tuple val(sample), path(bam)
    output: tuple val(sample), path(bam), path("${bam}.bai"), emit: bam
    script: "samtools index ${bam}"
}

// Optional UMI handling (enable with --umi; set --umi_pattern for your library).
process UMITOOLS_EXTRACT {
    label 'process_medium'
    container 'quay.io/biocontainers/umi_tools:1.1.5--py39hf95cd2a_0'
    input:  tuple val(sample), path(reads)
    output: tuple val(sample), path("${sample}.umi*.fastq.gz"), emit: reads
    script:
      if (reads instanceof List && reads.size() == 2)
        """
        umi_tools extract --bc-pattern=${params.umi_pattern} --bc-pattern2=${params.umi_pattern} \\
          -I ${reads[0]} --read2-in=${reads[1]} \\
          -S ${sample}.umi_1.fastq.gz --read2-out=${sample}.umi_2.fastq.gz
        """
      else
        """
        umi_tools extract --bc-pattern=${params.umi_pattern} -I ${reads} -S ${sample}.umi.fastq.gz
        """
}

process UMITOOLS_DEDUP {
    label 'process_medium'
    container 'quay.io/biocontainers/umi_tools:1.1.5--py39hf95cd2a_0'
    input:  tuple val(sample), path(bam), path(bai)
    output: tuple val(sample), path("${sample}.dedup.bam"), emit: bam
    script: "umi_tools dedup -I ${bam} -S ${sample}.dedup.bam"
}
