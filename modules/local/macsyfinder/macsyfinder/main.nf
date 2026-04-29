process MACSYFINDER_MACSYFINDER {
    tag "${meta.id}"
    label 'process_medium'

    container "oras://community.wave.seqera.io/library/macsyfinder:2.1.6--2680ae6a2ec885f7"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(faa)
    path(db)

    output:
    tuple val(meta), path("${prefix}_best_solution.tsv"), emit: tsv
    val(true), emit: done

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
    """
}