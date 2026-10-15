process HMMER_NHMMER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hdbdd923_1' :
        'biocontainers/hmmer:3.4--hdbdd923_1' }"

    input:
    tuple val(meta), path(fna)
    path(nod)
    path(tts)
    path(left)
    path(right)
    path(overdrive)

    output:
    tuple val(meta), path("${prefix}.nhmmer.all.tsv"), emit: tsv
    tuple val(meta), path("${prefix}.nodbox.tbl")  , emit: nodout
    tuple val(meta), path("${prefix}.ttsbox.tbl")  , emit: ttsout
    tuple val(meta), path("${prefix}.left.tbl")    , emit: leftout, optional: true
    tuple val(meta), path("${prefix}.right.tbl")   , emit: rightout, optional: true
    tuple val(meta), path("${prefix}.overdrive.tbl"), emit: odout, optional: true
    tuple val("${task.process}"), val('nhmmer'), eval("nhmmer -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//'"), emit: versions_nhmmer, topic: versions

    script:
    def args       = task.ext.args   ?: ''
    prefix         = task.ext.prefix ?: "${meta.id}"
    """
    nhmmer \\
        $args \\
        --tblout ${prefix}.nodbox.tbl \\
        --cpu $task.cpus \\
        $nod \\
        $fna
    
    nhmmer \\
        $args \\
        --tblout ${prefix}.ttsbox.tbl \\
        --cpu $task.cpus \\
        $tts \\
        $fna

    if [ -s "$left" ]; then
        nhmmer \\
            $args \\
            --tblout ${prefix}.left.tbl \\
            --cpu $task.cpus \\
            $left \\
            $fna

        nhmmer \\
            $args \\
            --tblout ${prefix}.right.tbl \\
            --cpu $task.cpus \\
            $right \\
            $fna

        nhmmer \\
            $args \\
            --tblout ${prefix}.overdrive.tbl \\
            --cpu $task.cpus \\
            $overdrive \\
            $fna
    fi

    echo "Parse and combine all nhmmer outputs"
    echo -e "contig\\tstart\\tend\\tstrand\\tevalue\\tscore\\ttype" > ${prefix}.nhmmer.all.tsv

    if [ -f "${prefix}.nodbox.tbl" ]; then
        grep -v "^#" ${prefix}.nodbox.tbl \\
        | awk '{print \$1"\\t"\$7"\\t"\$8"\\t"\$12"\\t"\$13"\\t"\$14"\\tnodbox"}' \\
        >> ${prefix}.nhmmer.all.tsv || true
    fi

    if [ -f "${prefix}.ttsbox.tbl" ]; then
        grep -v "^#" ${prefix}.ttsbox.tbl \\
        | awk '{print \$1"\\t"\$7"\\t"\$8"\\t"\$12"\\t"\$13"\\t"\$14"\\tttsbox"}' \\
        >> ${prefix}.nhmmer.all.tsv || true
    fi

    if [ -f "${prefix}.left.tbl" ]; then
        grep -v "^#" ${prefix}.left.tbl \\
        | awk '{print \$1"\\t"\$7"\\t"\$8"\\t"\$12"\\t"\$13"\\t"\$14"\\tleft"}' \\
        >> ${prefix}.nhmmer.all.tsv || true

        grep -v "^#" ${prefix}.right.tbl \\
        | awk '{print \$1"\\t"\$7"\\t"\$8"\\t"\$12"\\t"\$13"\\t"\$14"\\tright"}' \\
        >> ${prefix}.nhmmer.all.tsv || true

        grep -v "^#" ${prefix}.overdrive.tbl \\
        | awk '{print \$1"\\t"\$7"\\t"\$8"\\t"\$12"\\t"\$13"\\t"\$14"\\toverdrive"}' \\
        >> ${prefix}.nhmmer.all.tsv || true
    fi
    """

    stub:
    """
    touch ${prefix}.nhmmer.all.tsv
    touch ${prefix}.nodbox.tbl
    touch ${prefix}.ttsbox.tbl
    touch ${prefix}.left.tbl
    touch ${prefix}.right.tbl
    touch ${prefix}.overdrive.tbl
    """
}
