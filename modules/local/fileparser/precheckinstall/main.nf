process PRECHECK_INSTALL{
    publishDir "${params.beav_dir}", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d1/d1c10eb0e3654df8319ee582a6942f20a8b16c135b8662da53372f14ab57a79e/data'
        : 'community.wave.seqera.io/library/git:2.56.0--fe46a1319dd3d7e5'}"

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
    tuple val("${task.process}"), val('git'), eval("git -v 2>&1 | head -n 1"), emit: versions_git, topic: versions

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

    stub:
    """
    mkdir agrobacterium_taxa
    mkdir bakta_custom_protein
    mkdir databases
    mkdir models
    mkdir oncogenic_plasmids
    touch dif_db.fna
    touch more_T-DNA_borders.fna
    touch oriT_db.fna
    touch tiger.patch2
    """
}