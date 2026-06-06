#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { HARMONIZE }     from './modules/local/reference.nf'
include { FEATURECOUNTS } from './modules/local/quantify.nf'
include { FASTP; STAR_INDEX; STAR_ALIGN; SAMTOOLS_INDEX } from './modules/local/align.nf'
include { DESEQ2; PREDICT_G4; HYBRID_G4; G4_ENRICHMENT; GO_ENRICHMENT; GSEA; INTEGRATE } from './modules/local/analysis.nf'
include { REPORT }        from './modules/local/report.nf'

workflow {
    if (!params.input)     error "Provide --input <samplesheet.csv>"
    if (!params.fasta || !params.gff) error "Provide --fasta and --gff"
    if (!params.contrasts) error "Provide --contrasts <contrasts.csv>"

    def NO_FLAG   = file("${projectDir}/assets/NO_FLAG")
    def NO_MAP    = file("${projectDir}/assets/NO_MAP")
    def NO_BAM    = file("${projectDir}/assets/NO_BAM")
    def NO_G4     = file("${projectDir}/assets/NO_G4")
    def NO_HYBRID = file("${projectDir}/assets/NO_HYBRID")
    def NO_GSEA   = file("${projectDir}/assets/NO_GSEA")
    def header  = file(params.input).readLines()[0].split(',').collect{ it.trim() }
    def step    = params.step != 'auto' ? params.step :
                  ( params.count_matrix ? 'counts' : (header.contains('fastq_1') ? 'fastq' : 'bam') )
    log.info "nf-g4rnaseq | entry step = ${step}"

    samplesheet = file(params.input)
    contrasts   = file(params.contrasts)
    flag_genes  = params.flag_genes  ? file(params.flag_genes)  : NO_FLAG
    seqmap      = params.seqname_map ? file(params.seqname_map) : NO_MAP

    // BAM length-target for seqname harmonization (BAM entry only)
    def bam_target = NO_BAM
    if (step == 'bam' && params.harmonize == 'length') {
        def cols = file(params.input).readLines()[1].split(',')
        bam_target = file(cols[ header.indexOf('bam') ].trim())
    }

    // ---- reference: harmonize seqnames, build SAF + gene_info + indexed genome ----
    HARMONIZE(file(params.fasta), file(params.gff), bam_target, seqmap)
    genome    = HARMONIZE.out.fasta.first()
    gene_info = HARMONIZE.out.gene_info.first()
    saf       = HARMONIZE.out.saf.first()

    // ---- quantification ----
    if (step == 'counts') {
        if (!params.count_matrix) error "step=counts requires --count_matrix"
        counts = file(params.count_matrix)
    } else if (step == 'bam') {
        bam_ch  = Channel.fromPath(params.input).splitCsv(header:true).map { tuple(it.sample, file(it.bam)) }
        FEATURECOUNTS(bam_ch.map{it[1]}.collect(), bam_ch.map{it[0]}.collect(), saf)
        counts = FEATURECOUNTS.out.counts
    } else { // fastq
        reads = Channel.fromPath(params.input).splitCsv(header:true)
                  .map { r -> tuple(r.sample, r.fastq_2 ? [file(r.fastq_1), file(r.fastq_2)] : [file(r.fastq_1)]) }
        FASTP(reads)
        STAR_INDEX(genome)
        STAR_ALIGN(FASTP.out.reads, STAR_INDEX.out.index.first())
        SAMTOOLS_INDEX(STAR_ALIGN.out.bam)
        b = SAMTOOLS_INDEX.out.bam
        FEATURECOUNTS(b.map{it[1]}.collect(), b.map{it[0]}.collect(), saf)
        counts = FEATURECOUNTS.out.counts
    }

    // ---- differential expression ----
    DESEQ2(counts, samplesheet, gene_info, contrasts, flag_genes)
    deseq2 = DESEQ2.out.dir.first()

    // ---- G4 + hybrid-G4 (optional) ----
    g4_dir     = Channel.value(NO_G4)
    hybrid_dir = Channel.value(NO_HYBRID)
    if (!params.skip_g4)     { PREDICT_G4(genome, gene_info); g4_dir     = PREDICT_G4.out.dir.first() }
    if (!params.skip_hybrid) { HYBRID_G4(genome, gene_info);  hybrid_dir = HYBRID_G4.out.dir.first() }

    // ---- G4<->DE enrichment / GO / GSEA / integrate ----
    if (!params.skip_g4 && !params.skip_hybrid) G4_ENRICHMENT(deseq2, g4_dir, hybrid_dir, flag_genes)
    GO_ENRICHMENT(deseq2, flag_genes)
    gsea_dir = Channel.value(NO_GSEA)
    if (!params.skip_gsea && !params.skip_g4 && !params.skip_hybrid) {
        GSEA(deseq2, g4_dir, hybrid_dir); gsea_dir = GSEA.out.dir.ifEmpty(NO_GSEA).first()
    }
    INTEGRATE(gene_info, g4_dir, hybrid_dir, deseq2, flag_genes)

    if (!params.skip_report) REPORT(deseq2, g4_dir, hybrid_dir, gsea_dir, INTEGRATE.out.dir.first())
}
