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
    tuple val("${task.process}"), val('defensefinder'), eval("defense-finder --version 2>&1 | head -n 1"), topic: versions, emit: versions_defensefinder

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    export HOME=\$PWD/.defensefinder_home
    mkdir -p "\$HOME"

    export OUTDIR=\$PWD/output

    defense-finder run --db-type ordered_replicon \\
                       --out-dir \$OUTDIR \\
                       --models-dir ${db} \\
                       ${faa}
                       
    mv \$OUTDIR/${prefix}_defense_finder_genes.tsv .
    """

    stub:
    """
    touch ${prefix}_defense_finder_genes.tsv
    """
}