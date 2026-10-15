process DL_MACSY_DB{
    label 'process_medium'

    cache true
    publishDir "${params.beav_dir}/models", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/28/28cd67a8ad35c2b4cbade5bc2e5533ab7ad3dc39843d5c6e8dd54d84518055c2/data'
        : 'community.wave.seqera.io/library/macsyfinder_biopython_pandas_python_requests:1cfad224e5c82c97'}"

    output:
    path("macsy-finder-models"), emit: db
    tuple val("${task.process}"), val('macsyfinder'), eval("macsyfinder --version 2>&1 | head -n 1"), emit: versions_macsyfinder, topic: versions
    tuple val("${task.process}"), val('macsydata'), eval("macsydata --version 2>&1 | head -n 1"), emit: versions_macsydata, topic: versions

    script:
    """
    macsydata install --target macsy-finder-models -U TXSScan
    """

    stub:
    """
    mkdir macsy-finder-models
    """
}