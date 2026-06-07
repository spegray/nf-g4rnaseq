#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { HARMONIZE }     from './modules/local/reference.nf'
include { FEATURECOUNTS } from './modules/local/quantify.nf'
include { FASTP; STAR_INDEX; STAR_ALIGN; SAMTOOLS_INDEX; UMITOOLS_EXTRACT; UMITOOLS_DEDUP } from './modules/local/align.nf'
include { DESEQ2; PREDICT_G4; HYBRID_G4; G4_ENRICHMENT; GO_ENRICHMENT; GSEA; INTEGRATE } from './modules/local/analysis.nf'
include { REPORT }        from './modules/local/report.nf'
include { MULTIQC; VERSIONS } from './modules/local/qc.nf'

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
    def NO_G4ENRICH = file("${projectDir}/assets/NO_G4ENRICH")
    def NO_GO     = file("${projectDir}/assets/NO_GO")
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
        // Pair [sample, bam] and sort ONCE, then derive the bam list and the name
        // list from that same ordered structure, so featureCounts columns are
        // always labeled with the correct sample (no independent-collect order drift).
        ordered = Channel.fromPath(params.input).splitCsv(header:true)
                    .map { tuple(it.sample, file(it.bam)) }
                    .toSortedList { a, b -> a[0] <=> b[0] }
        FEATURECOUNTS(ordered.map{ it.collect{ t -> t[1] } }, ordered.map{ it.collect{ t -> t[0] } }, saf)
        counts = FEATURECOUNTS.out.counts
    } else { // fastq
        reads = Channel.fromPath(params.input).splitCsv(header:true)
                  .map { r -> tuple(r.sample, r.fastq_2 ? [file(r.fastq_1), file(r.fastq_2)] : [file(r.fastq_1)]) }
        FASTP(reads)
        STAR_INDEX(genome)
        if (params.umi) {                       // optional UMI extract -> align -> dedup
            UMITOOLS_EXTRACT(FASTP.out.reads)
            STAR_ALIGN(UMITOOLS_EXTRACT.out.reads, STAR_INDEX.out.index.first())
            SAMTOOLS_INDEX(STAR_ALIGN.out.bam)
            UMITOOLS_DEDUP(SAMTOOLS_INDEX.out.bam)
            aln = UMITOOLS_DEDUP.out.bam
        } else {
            STAR_ALIGN(FASTP.out.reads, STAR_INDEX.out.index.first())
            SAMTOOLS_INDEX(STAR_ALIGN.out.bam)
            aln = SAMTOOLS_INDEX.out.bam.map { t -> tuple(t[0], t[1]) }   // drop the .bai
        }
        ordered = aln.toSortedList { a, b -> a[0] <=> b[0] }
        FEATURECOUNTS(ordered.map{ it.collect{ t -> t[1] } }, ordered.map{ it.collect{ t -> t[0] } }, saf)
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
    g4enrich_dir = Channel.value(NO_G4ENRICH)
    if (!params.skip_g4 && !params.skip_hybrid) {
        G4_ENRICHMENT(deseq2, g4_dir, hybrid_dir, flag_genes)
        g4enrich_dir = G4_ENRICHMENT.out.dir.ifEmpty(NO_G4ENRICH).first()
    }
    GO_ENRICHMENT(deseq2, flag_genes)
    go_dir = GO_ENRICHMENT.out.dir.ifEmpty(NO_GO).first()
    gsea_dir = Channel.value(NO_GSEA)
    if (!params.skip_gsea && !params.skip_g4 && !params.skip_hybrid) {
        GSEA(deseq2, g4_dir, hybrid_dir); gsea_dir = GSEA.out.dir.ifEmpty(NO_GSEA).first()
    }
    INTEGRATE(gene_info, g4_dir, hybrid_dir, deseq2, flag_genes)

    if (!params.skip_report)
        REPORT(deseq2, g4_dir, hybrid_dir, g4enrich_dir, go_dir, gsea_dir, INTEGRATE.out.dir.first(),
               file("${projectDir}/assets/report_template.qmd"))

    // ---- provenance + read/alignment QC aggregation ----
    VERSIONS()
    if (!params.skip_multiqc) {
        if (step == 'fastq')
            MULTIQC( FASTP.out.json.mix(STAR_ALIGN.out.log).mix(FEATURECOUNTS.out.summary).collect() )
        else if (step == 'bam')
            MULTIQC( FEATURECOUNTS.out.summary.collect() )
    }
}
