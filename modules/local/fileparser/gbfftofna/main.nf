process CONVERT_GBFF_TO_FNA{

    input:
    tuple val(meta), path(gbff)
    val(done)

    output:
    path("${meta.id}.fna"), emit: fna
    val(true), emit: done

    script:
    """
python - <<'PYCODE'

from Bio import SeqIO

SeqIO.convert("${gbff}", "genbank", "${gbff}.fna", "fasta")

PYCODE
    """

}