process GSEA_HEATMAP {
    label 'process_low'

    conda "bioconda::bioconductor-deseq2=1.34.0"
    container "biocontainers-r-custom-ComplexHeatmap"

    input:
    path(input_files)

    output:
    path "*heatmap.png"                                        , emit: heatmap
    path "R_sessionInfo.log"                                   , emit: session_info
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'heatmap.R'
}