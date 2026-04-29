process PRECHECK_INSTALL{
    publishDir "${params.beav_dir}", mode: 'copy'

    input:
    path(local)

    output:
    path("agrobacterium_taxa/")
    path("bakta_custom_protein/")
    path("databases/")
    path("models/")
    path("oncogenic_plasmids/")
    path("dif_db.fna")
    path("more_T-DNA_borders.fna")
    path("oriT_db.fna")
    path("tiger.patch2")
    val(true), emit: done

    script:
    """
    set -e

    GITHUB_REPO="https://github.com/weisberglab/beavDB.git"

    rm -rf ${local}

    git clone --depth 1 \$GITHUB_REPO ${local}

    mkdir -p ${local}/databases
    mkdir -p ${local}/models/TXSS
    mv ${local}/* .
    """
}