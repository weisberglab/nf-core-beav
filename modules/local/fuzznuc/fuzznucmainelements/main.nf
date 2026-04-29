process FUZZNUC_MAIN_ELEMENTS{
    tag "${meta.id}"
    label 'process_medium'

    container "oras://community.wave.seqera.io/library/emboss:6.6.0--b0e4491212918760"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fna)
    val(agro)

    output:
    tuple val(meta), path("${prefix}.combined.tsv"), emit: tsv
    tuple val(meta), path("${prefix}.hrp_box.out") , emit: hrp
    tuple val(meta), path("${prefix}.pip_box.out") , emit: pip
    tuple val(meta), path("${prefix}.kops.out")    , emit: kops
    tuple val(meta), path("${prefix}.matS.out")    , emit: matS
    tuple val(meta), path("${prefix}.virbox.out")  , emit: vir, optional: true
    tuple val(meta), path("${prefix}.trabox.out")  , emit: tra, optional: true
    val(true), emit: done

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    fuzznuc -sequence ${fna} \\
            -pattern 'GGAAC[CT]N(15,17)CCACNNA' \\
            -complement \\
            -rformat excel \\
            -outfile ${prefix}.hrp_box.out


    fuzznuc -sequence ${fna} \\
            -pattern 'TTCGBN(15)TTCGB' \\
            -complement \\
            -rformat excel \\
            -outfile ${prefix}.pip_box.out

    fuzznuc -sequence ${fna} \\
            -pattern 'GGGNAGGG' \\
            -complement \\
            -rformat excel \\
            -outfile ${prefix}.kops.out

    fuzznuc -sequence ${fna} \\
            -pattern 'GTGACANTGTCAC' \\
            -complement \\
            -rformat excel \\
            -outfile ${prefix}.matS.out

    if [ "${agro}" == "true" ]; then
            fuzznuc -sequence ${fna} \\
                -pattern 'RTTDCAWWTGHAAY' \\
                -rformat excel \\
                -complement \\
                -outfile ${prefix}.virbox.out
    
        fuzznuc -sequence ${fna} \\
                -pattern 'WNGTGMARAWYTGCACDW' \\
                -rformat excel \\
                -complement \\
                -outfile ${prefix}.trabox.out
    fi
        

    # Combine outputs into one TSV with Motif column
    awk 'NR==1{print \$0"\\tMotif"} NR>1{print \$0"\\tHRP"}' ${prefix}.hrp_box.out > ${prefix}.hrp_tagged.tsv
    awk 'NR==1{print \$0"\\tMotif"} NR>1{print \$0"\\tPIP"}' ${prefix}.pip_box.out > ${prefix}.pip_tagged.tsv
    awk 'NR==1{print \$0"\\tMotif"} NR>1{print \$0"\\tKOPS"}' ${prefix}.kops.out > ${prefix}.kops_tagged.tsv
    awk 'NR==1{print \$0"\\tMotif"} NR>1{print \$0"\\tmatS"}' ${prefix}.matS.out > ${prefix}.matS_tagged.tsv

    if [ "${agro}" == "true" ]; then
        awk 'NR==1{print \$0"\\tMotif"} NR>1{print \$0"\\tVIR"}' ${prefix}.virbox.out > ${prefix}.vir_tagged.tsv
        awk 'NR==1{print \$0"\\tMotif"} NR>1{print \$0"\\tTRA"}' ${prefix}.trabox.out > ${prefix}.tra_tagged.tsv

        # Merge all non-empty tagged files
        cat ${prefix}.hrp_tagged.tsv ${prefix}.pip_tagged.tsv ${prefix}.kops_tagged.tsv ${prefix}.matS_tagged.tsv ${prefix}.vir_tagged.tsv ${prefix}.tra_tagged.tsv \\
                | awk 'NR==1 || FNR>1' \\
                | sort -k1,1 -k2,2n > ${prefix}.combined.tsv
    else
        # Merge all non-empty tagged files
        cat ${prefix}.hrp_tagged.tsv ${prefix}.pip_tagged.tsv ${prefix}.kops_tagged.tsv ${prefix}.matS_tagged.tsv \\
                | awk 'NR==1 || FNR>1' \\
                | sort -k1,1 -k2,2n > ${prefix}.combined.tsv
    fi


    """

}