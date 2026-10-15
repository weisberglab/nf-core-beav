process DL_TXSS_PAPBLAST{
    errorStrategy 'retry'
    maxRetries 3
    label 'process_medium'

    publishDir "${params.beav_dir}/models", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1d/1d9f21eaec6d6273b80db8ab82bba360ddfe6fdc647d44814032323f671496bd/data'
        : 'community.wave.seqera.io/library/defense-finder_python:9760d1f62683751b'}"

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