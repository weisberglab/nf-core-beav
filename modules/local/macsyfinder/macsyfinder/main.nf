process MACSYFINDER_MACSYFINDER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/28/28cd67a8ad35c2b4cbade5bc2e5533ab7ad3dc39843d5c6e8dd54d84518055c2/data'
        : 'community.wave.seqera.io/library/macsyfinder_biopython_pandas_python_requests:1cfad224e5c82c97'}"

    input:
    tuple val(meta), path(faa)
    path(db)

    output:
    tuple val(meta), path("${prefix}_best_solution.tsv"), emit: tsv
    path "versions.yml"                             , emit: versions

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p results

    macsyfinder --db-type ordered_replicon \\
                --sequence-db ${faa} \\
                --models-dir ${db} \\
                --models TXSScan \\
                --out-dir ./results

    mv ./results/best_solution.tsv ${prefix}_best_solution.tsv

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        macsyfinder: \$(macsyfinder --version 2>&1 | head -n 1)
        macsydata: \$(macsydata --version 2>&1 | head -n 1)
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}_best_solution.tsv

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        macsyfinder: \$(macsyfinder --version 2>&1 | head -n 1)
        macsydata: \$(macsydata --version 2>&1 | head -n 1)
    END_VERSIONS
    """
}