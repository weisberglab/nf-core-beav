process BLAST_BLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0c/0c86cbb145786bf5c24ea7fb13448da5f7d5cd124fd4403c1da5bc8fc60c2588/data':
        'oras://community.wave.seqera.io/library/bedtools_blast:e4f7c8c6d31f67f9' }"
    // ""
    
    input:
    tuple val(meta) , path(fasta)
    path oriT
    path dif
    path taxidlist
    val taxids
    val negative_tax

    output:
    tuple val(meta), path("${prefix}_oriT.txt"), path("${prefix}_difSites.txt"), emit: out
    path "versions.yml"                                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // def is_compressed = fasta.getExtension() == "gz" ? true : false
    // def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    def negative_tax_cmd = negative_tax ? "negative_" : ""
    def taxidlist_cmd = taxidlist ? "-${negative_tax_cmd}taxidlist ${taxidlist}" : ""
    def taxids_cmd = taxids ? "-${negative_tax_cmd}taxids ${taxids}" : ""
    if (taxidlist_cmd.any() && taxids_cmd.any()) {
        log.error("ERROR: taxidlist and taxids can not be used at the same time, choose only one argument to use for tax id filtering.")
    }

    """
    makeblastdb -in ${fasta} -dbtype nucl 1> /dev/null

    blastn \\
        -task blastn-short \\
        -outfmt '6 std qlen slen qseq sseq' \\
        -num_threads ${task.cpus} \\
        -dust no \\
        -db ${fasta} \\
        -query ${oriT} \\
        ${taxidlist_cmd} \\
        ${taxids_cmd} \\
        ${args} \\
        -out ${prefix}_oriT.txt \\

    blastn \\
        -task blastn-short \\
        -outfmt '6 std qlen slen qseq sseq' \\
        -num_threads ${task.cpus} \\
        -dust no \\
        -db ${fasta} \\
        -query ${dif} \\
        ${taxidlist_cmd} \\
        ${taxids_cmd} \\
        ${args} \\
        -out ${prefix}_difSites.txt \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}_oriT.txt
    touch ${prefix}_difSites.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}