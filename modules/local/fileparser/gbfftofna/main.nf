process CONVERT_GBFF_TO_FNA{

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/20/20ae2f5b7f02eccb8e16aec73e13ba03b11744576ae70698e5028bc806193c9c/data'
        : 'community.wave.seqera.io/library/biopython_seqiolib:63d1dbe8c336c0d9'}"

    input:
    tuple val(meta), path(gbff)

    output:
    path("${meta.id}.fna"), emit: fna
    tuple val("${task.process}"), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), emit: versions_biopython, topic: versions

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
python - <<'PYCODE'

from Bio import SeqIO

SeqIO.convert("${gbff}", "genbank", "${gbff}.fna", "fasta")

PYCODE
    """

    stub:
    """
    touch ${prefix}.fna
    """

}