process DEFENSEFINDER_DEFENSEFINDER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1d/1d9f21eaec6d6273b80db8ab82bba360ddfe6fdc647d44814032323f671496bd/data'
        : 'community.wave.seqera.io/library/defense-finder_python:9760d1f62683751b'}"

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