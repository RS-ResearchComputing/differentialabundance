process GEOQUERY_GETGEO {
    tag "$meta.id"
    label 'process_single'

    conda 'modules/nf-core/geoquery/getgeo/environment.yml'
    container "geoquery:2.66.0--r42hdfd78af_0"

    input:
    tuple val(meta), val(querygse)

    output:
    tuple val(meta), path("*.rds")            , emit: rds
    tuple val(meta), path("*matrix.tsv")      , emit: expression
    tuple val(meta), path("*annotation.tsv")  , emit: annotation
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'getgeo.R'
}
