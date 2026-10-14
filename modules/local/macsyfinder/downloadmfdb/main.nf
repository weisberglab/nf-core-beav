process DL_MACSY_DB{
    label 'process_medium'

    cache true
    publishDir "${params.beav_dir}/models", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/macsyfinder:2.1.6--2680ae6a2ec885f7"

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