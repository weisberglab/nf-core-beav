process DEFENSEFINDER_DEFENSEFINDER {
    tag "${meta.id}"
    label 'process_medium'

    container "quay.io/biocontainers/defense-finder:2.0.1--pyhdfd78af_0"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(faa)
    path(db)

    output:
    tuple val(meta), path("${prefix}_defense_finder_genes.tsv"), emit: tsv
    val(true), emit: done

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    export HOME=\$PWD/.defensefinder_home
    mkdir -p "\$HOME"

    defense-finder run --db-type ordered_replicon \\
                       --out-dir . \\
                       --models-dir ${db} \\
                       ${faa}
    """
}