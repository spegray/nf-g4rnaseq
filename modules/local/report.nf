// Assemble a summary report from the pipeline results.
process REPORT {
    label 'process_low'
    container params.r_container
    input:
      path deseq2
      path g4         // NO_FILE if --skip_g4
      path hybrid     // NO_FILE if --skip_hybrid
      path gsea       // NO_FILE if skipped
      path integrate
    output:
      path "summary_report.md", emit: report
    script:
      def g4o  = g4.name     != 'NO_G4'     ? "--g4_dir ${g4}"          : ''
      def hyo  = hybrid.name != 'NO_HYBRID' ? "--hybrid_dir ${hybrid}"  : ''
      def gso  = gsea.name   != 'NO_GSEA'   ? "--gsea_dir ${gsea}"      : ''
      """
      render_report.R --deseq2_dir ${deseq2} ${g4o} ${hyo} ${gso} --integrate_dir ${integrate} --outdir .
      """
}
