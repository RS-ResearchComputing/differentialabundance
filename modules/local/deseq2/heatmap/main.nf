process DESEQ2_HEATMAP {
    tag "$meta"
    label 'process_low'

    conda "bioconda::bioconductor-deseq2=1.34.0"
    container "biocontainers-r-custom-ComplexHeatmap"

    input:
    path(input_files)
    tuple val(logFC_column), val(logFC_threshold)
    tuple val(padj_column), val(padj_threshold)
    tuple val(meta), path(contrasts)

    output:
    path "*heatmap.png"                                        , emit: heatmap
    path "R_sessionInfo.log"                                 , emit: session_info
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'heatmap.R'
}