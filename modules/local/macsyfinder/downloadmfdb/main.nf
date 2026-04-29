process DL_MACSY_DB{
    label 'process_medium'

    cache true
    publishDir "${params.beav_dir}/models", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/macsyfinder:2.1.6--2680ae6a2ec885f7"

    input:
    val(dummy)

    output:
    path("macsy-finder-models"), emit: db

    script:
    """
    macsydata install --target macsy-finder-models -U TXSScan
    """
}