process CONVERT_GBFF_TO_FNA{

    input:
    tuple val(meta), path(gbff)
    val(done)

    output:
    path("${meta.id}.fna"), emit: fna

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