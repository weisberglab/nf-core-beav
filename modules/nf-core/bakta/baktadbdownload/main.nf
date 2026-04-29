process BAKTA_BAKTADBDOWNLOAD {
    tag "${params.bakta_db_type}"
    label 'process_single'
    
    publishDir "${params.beav_dir}/databases", mode: 'move', pattern: '*db'

    conda "${moduleDir}/environment.yml"
    container "https://depot.galaxyproject.org/singularity/bakta:1.12.0--pyhdfd78af_0"

    input:
    val(dummy)
    
    output:
    path "*db", emit: db
    path "versions.yml", emit: versions
    val(true),  emit: done

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bakta_db \\
        download \\
        ${args} \\
        --type ${params.bakta_db_type} || true

    if [[ ${params.bakta_db_type} == 'light' ]]; then
        mv db-light bakta_db
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta_db --version) 2>&1 | cut -f '2' -d ' ')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo "bakta_db \\
        download \\
        ${args}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta_db --version) 2>&1 | cut -f '2' -d ' ')
    END_VERSIONS
    """
}
