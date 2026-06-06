// Render the parameterized Quarto report (HTML) from the pipeline results.
// Falls back to a lightweight Markdown summary if Quarto is unavailable or errors.
process REPORT {
    label 'process_low'
    container params.r_container
    input:
      path deseq2
      path g4
      path hybrid
      path g4enrich
      path go
      path gsea
      path integrate
      path qmd
    output:
      path "report.html",       emit: html, optional: true
      path "summary_report.md", emit: md,   optional: true
    script:
      def g4o = g4.name     != 'NO_G4'     ? "--g4_dir ${g4}"         : ''
      def hyo = hybrid.name != 'NO_HYBRID' ? "--hybrid_dir ${hybrid}" : ''
      def gso = gsea.name   != 'NO_GSEA'   ? "--gsea_dir ${gsea}"     : ''
      """
      # Quarto reads the staged result dirs (deseq2/, g4/, hybrid/, g4enrich/,
      # go/, gsea/, integrate/) by name; missing/placeholder dirs are skipped.
      cp ${qmd} report.qmd
      quarto render report.qmd --to html > quarto_render.log 2>&1 \\
        || echo "[REPORT] quarto render failed (see quarto_render.log) -> markdown summary only"
      # always also emit the dependency-free markdown summary (and as a fallback)
      render_report.R --deseq2_dir ${deseq2} ${g4o} ${hyo} ${gso} --integrate_dir ${integrate} --outdir . || true
      """
}
