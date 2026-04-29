process DL_TXSS_PAPBLAST{
    errorStrategy 'retry'
    maxRetries 5
    label 'process_medium'

    publishDir "${params.beav_dir}/models", mode: 'copy'

    container "quay.io/biocontainers/defense-finder:2.0.1--pyhdfd78af_0"    
    conda "${moduleDir}/environment.yml"

    input:
    val(dummy)
    
    output:
    path("defense-finder-models"), emit: db

    script:
    """
    defense-finder update --models-dir "defense-finder-models"

    # defense finder hates me because it installs a version of CasFinder that is too new for the macsyfinder version it uses
    cd defense-finder-models
    rm -rf CasFinder
    macsydata install CasFinder==3.1.0 --models-dir tmpc
    mv tmpc/CasFinder .
    rm -rf tmpc
    cd ..
    """
}