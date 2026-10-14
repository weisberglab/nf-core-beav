process DL_TXSS_PAPBLAST{
    errorStrategy 'retry'
    maxRetries 5
    label 'process_medium'

    publishDir "${params.beav_dir}/models", mode: 'copy'

    container "quay.io/biocontainers/defense-finder:2.0.1--pyhdfd78af_0"    
    conda "${moduleDir}/environment.yml"
    
    output:
    path("defense-finder-models"), emit: db
    tuple val("${task.process}"), val('defense-finder'), eval("defense-finder --version 2>&1 | head -n 1"), emit: versions_defensefinder, topic: versions

    script:
    """
    defense-finder update --models-dir "defense-finder-models"

    cd defense-finder-models
    rm -rf CasFinder
    macsydata install CasFinder==3.1.0 --models-dir tmpc
    mv tmpc/CasFinder .
    rm -rf tmpc
    cd ..
    """

    stub:
    """
    mkdir defense-finder-models
    """
}