process HMMER_NHMMER {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hdbdd923_1' :
        'biocontainers/hmmer:3.4--hdbdd923_1' }"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fna)
    path(nod)
    path(tts)
    path(left)
    path(right)
    path(overdrive)

    output:
    tuple val(meta), path("${meta.id}.nhmmer.all.tsv"), emit: tsv
    tuple val(meta), path("${meta.id}.nodbox.tbl")  , emit: nodout
    tuple val(meta), path("${meta.id}.ttsbox.tbl")  , emit: ttsout
    tuple val(meta), path("${meta.id}.left.tbl")    , emit: leftout, optional: true
    tuple val(meta), path("${meta.id}.right.tbl")   , emit: rightout, optional: true
    tuple val(meta), path("${meta.id}.overdrive.tbl"), emit: odout, optional: true
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ''
    prefix         = task.ext.prefix ?: "${meta.id}"
    """
    nhmmer \\
        $args \\
        --tblout ${meta.id}.nodbox.tbl \\
        --cpu $task.cpus \\
        $nod \\
        $fna
    
    nhmmer \\
        $args \\
        --tblout ${meta.id}.ttsbox.tbl \\
        --cpu $task.cpus \\
        $tts \\
        $fna

    if [ -s "$left" ]; then
        nhmmer \\
            $args \\
            --tblout ${meta.id}.left.tbl \\
            --cpu $task.cpus \\
            $left \\
            $fna

        nhmmer \\
            $args \\
            --tblout ${meta.id}.right.tbl \\
            --cpu $task.cpus \\
            $right \\
            $fna

        nhmmer \\
            $args \\
            --tblout ${meta.id}.overdrive.tbl \\
            --cpu $task.cpus \\
            $overdrive \\
            $fna
    fi

    echo "Parse and combine all nhmmer outputs"
    echo -e "contig\\tstart\\tend\\tstrand\\tevalue\\tscore\\ttype" > ${meta.id}.nhmmer.all.tsv

    if [ -f "${meta.id}.nodbox.tbl" ]; then
        grep -v "^#" ${meta.id}.nodbox.tbl \\
        | awk '{print \$1"\\t"\$7"\\t"\$8"\\t"\$12"\\t"\$13"\\t"\$14"\\tnodbox"}' \\
        >> ${meta.id}.nhmmer.all.tsv || true
    fi

    if [ -f "${meta.id}.ttsbox.tbl" ]; then
        grep -v "^#" ${meta.id}.ttsbox.tbl \\
        | awk '{print \$1"\\t"\$7"\\t"\$8"\\t"\$12"\\t"\$13"\\t"\$14"\\tttsbox"}' \\
        >> ${meta.id}.nhmmer.all.tsv || true
    fi

    if [ -f "${meta.id}.left.tbl" ]; then
        grep -v "^#" ${meta.id}.left.tbl \\
        | awk '{print \$1"\\t"\$7"\\t"\$8"\\t"\$12"\\t"\$13"\\t"\$14"\\tleft"}' \\
        >> ${meta.id}.nhmmer.all.tsv || true

        grep -v "^#" ${meta.id}.right.tbl \\
        | awk '{print \$1"\\t"\$7"\\t"\$8"\\t"\$12"\\t"\$13"\\t"\$14"\\tright"}' \\
        >> ${meta.id}.nhmmer.all.tsv || true

        grep -v "^#" ${meta.id}.overdrive.tbl \\
        | awk '{print \$1"\\t"\$7"\\t"\$8"\\t"\$12"\\t"\$13"\\t"\$14"\\toverdrive"}' \\
        >> ${meta.id}.nhmmer.all.tsv || true
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(nhmmer -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """
}
