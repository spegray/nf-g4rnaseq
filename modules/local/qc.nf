// Read/alignment QC aggregation (MultiQC) and software-version provenance.

process MULTIQC {
    label 'process_low'
    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'
    input:
      path 'logs/*'        // featureCounts summaries, STAR Log.final, fastp json, ...
    output:
      path "multiqc_report.html", emit: report, optional: true
      path "multiqc_data",        emit: data,   optional: true
    script:
      """
      multiqc logs --filename multiqc_report.html --force || echo "[MULTIQC] nothing to aggregate"
      """
}

process VERSIONS {
    label 'process_low'
    container params.r_container
    output:
      path "software_versions.yml"
    script:
      """
      {
        echo "nf-g4rnaseq: ${workflow.manifest.version}"
        echo "nextflow: ${workflow.nextflow.version}"
        command -v samtools      >/dev/null 2>&1 && echo "samtools: \$(samtools --version 2>/dev/null | head -1 | awk '{print \$NF}')"
        command -v featureCounts >/dev/null 2>&1 && echo "subread: \$(featureCounts -v 2>&1 | grep -oE 'v[0-9.]+')"
        command -v STAR          >/dev/null 2>&1 && echo "STAR: \$(STAR --version 2>/dev/null)"
        command -v fastp         >/dev/null 2>&1 && echo "fastp: \$(fastp --version 2>&1 | awk '{print \$NF}')"
        Rscript -e 'cat(sprintf("R: %s\\n", getRversion())); for (p in c("DESeq2","apeglm","ashr","pqsfinder","Biostrings","GenomicRanges","rtracklayer","clusterProfiler","limma","fgsea")) if (requireNamespace(p, quietly=TRUE)) cat(sprintf("%s: %s\\n", p, as.character(packageVersion(p))))' 2>/dev/null
      } > software_versions.yml
      """
}
