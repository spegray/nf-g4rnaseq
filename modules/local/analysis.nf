// ===========================================================================
// Core analysis processes (wrap the validated bin/ R scripts). Each writes to a
// named output directory and emits it, keeping channel wiring simple.
// The R container is used under -profile docker/singularity; under the local
// profile the bin/ scripts (auto-added to PATH by Nextflow) run with PATH tools.
// ===========================================================================

process DESEQ2 {
    label 'process_medium'
    container params.r_container
    input:
      path counts
      path samplesheet
      path gene_info
      path contrasts
      path flag_genes        // assets/NO_FILE if none
    output:
      path "deseq2", emit: dir
    script:
      def fg = flag_genes.name != 'NO_FLAG' ? "--flag_genes ${flag_genes}" : ''
      def gi = gene_info.name != 'NO_FILE'  ? "--gene_info ${gene_info}"   : ''
      """
      mkdir -p deseq2
      run_deseq2.R --counts ${counts} --samplesheet ${samplesheet} ${gi} \\
        --design '${params.design_formula}' --contrasts ${contrasts} \\
        --reference_levels '${params.reference_levels}' ${fg} \\
        --min_count ${params.min_count} --padj ${params.padj} --outdir deseq2
      """
}

process PREDICT_G4 {
    label 'process_medium'
    container params.r_container
    input:
      path fasta
      path gene_info
    output:
      path "g4", emit: dir
    script:
      """
      mkdir -p g4
      predict_g4.R --fasta ${fasta} --gene_info ${gene_info} \\
        --min_score ${params.g4_min_score} --highconf_score ${params.g4_highconf_score} \\
        --promoter ${params.g4_promoter} --outdir g4
      """
}

process HYBRID_G4 {
    label 'process_medium'
    container params.r_container
    input:
      path fasta
      path gene_info
    output:
      path "hybrid", emit: dir
    script:
      """
      mkdir -p hybrid
      hybrid_g4.R --fasta ${fasta} --gene_info ${gene_info} \\
        --run_len ${params.hybrid_run_len} --loop_max ${params.hybrid_loop_max} \\
        --span_max ${params.hybrid_span_max} --min_tracts ${params.hybrid_min_tracts} \\
        --min_expr ${params.min_expr} --outdir hybrid
      """
}

process G4_ENRICHMENT {
    label 'process_low'
    container params.r_container
    input:
      path deseq2
      path g4
      path hybrid
      path flag_genes
    output:
      path "g4enrich", emit: dir
    script:
      def fg = flag_genes.name != 'NO_FLAG' ? "--flag_genes ${flag_genes}" : ''
      def od = params.orgdb ? "--orgdb ${params.orgdb} --keytype ${params.gene_id_keytype}" : ''
      """
      mkdir -p g4enrich
      g4_enrichment.R --deseq2_dir ${deseq2} --g4_features ${g4}/gene_g4_features.tsv \\
        --hybrid ${hybrid}/gene_hybrid_g4.tsv ${fg} ${od} --padj ${params.padj} --outdir g4enrich
      """
}

process GO_ENRICHMENT {
    label 'process_low'
    container params.r_container
    input:
      path deseq2
      path flag_genes
    output:
      path "go", emit: dir, optional: true
    when:
      !params.skip_go && params.orgdb != null
    script:
      def fg = flag_genes.name != 'NO_FLAG' ? "--flag_genes ${flag_genes}" : ''
      def kg = params.kegg_organism ? "--kegg_organism ${params.kegg_organism}" : ''
      """
      mkdir -p go
      go_enrichment.R --deseq2_dir ${deseq2} --orgdb ${params.orgdb} --keytype ${params.gene_id_keytype} \\
        ${kg} ${fg} --padj ${params.padj} --outdir go
      """
}

process GSEA {
    label 'process_low'
    container params.r_container
    input:
      path deseq2
      path g4
      path hybrid
    output:
      path "gsea", emit: dir, optional: true
    when:
      !params.skip_gsea
    script:
      def od = params.orgdb ? "--orgdb ${params.orgdb} --keytype ${params.gene_id_keytype}" : ''
      """
      mkdir -p gsea
      gsea.R --deseq2_dir ${deseq2} --g4_features ${g4}/gene_g4_features.tsv \\
        --hybrid ${hybrid}/gene_hybrid_g4.tsv ${od} --outdir gsea
      """
}

process INTEGRATE {
    label 'process_low'
    container params.r_container
    input:
      path gene_info
      path g4
      path hybrid
      path deseq2
      path flag_genes
    output:
      path "integrate", emit: dir
    script:
      def fg  = flag_genes.name != 'NO_FLAG'   ? "--flag_genes ${flag_genes}"               : ''
      def g4o = g4.name        != 'NO_G4'      ? "--g4_features ${g4}/gene_g4_features.tsv" : ''
      def hyo = hybrid.name    != 'NO_HYBRID'  ? "--hybrid ${hybrid}/gene_hybrid_g4.tsv"    : ''
      """
      mkdir -p integrate
      integrate.R --gene_info ${gene_info} ${g4o} ${hyo} --deseq2_dir ${deseq2} ${fg} \\
        --padj ${params.padj} --outdir integrate
      """
}
