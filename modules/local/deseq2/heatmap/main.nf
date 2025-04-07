process DESEQ2_HEATMAP {
    tag "$meta2"
    label 'process_medium'

    conda "bioconda::bioconductor-deseq2=1.34.0"
    container "biocontainers-r-GRCm39-GRCh38"

    input:
    path(input_files)
    tuple val(logFC_column), val(logFC_threshold)
    tuple val(padj_column), val(padj_threshold)
    tuple val(meta), path(contrasts)
    val(genome)
    tuple val(meta2), path(gtf_file)

    output:
    path "*heatmap.png"                                        , emit: heatmap
    path "R_sessionInfo.log"                                 , emit: session_info
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'heatmap.R'
}