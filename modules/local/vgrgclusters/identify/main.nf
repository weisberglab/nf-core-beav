process IDENTIFY_VGRG_CLUSTERS{
    label 'process_medium'
    tag "$meta.id"

    container "oras://community.wave.seqera.io/library/hmmer:3.4--776668bdf0589c68"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path('*')
    path(hmm1)
    path(hmm2)

    output:
    tuple val(meta), path("${prefix}_T6SS_vgrG_cluster_list.out"), emit: clist

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    {
        if [[ -f "${prefix}.gff3" ]]; then
            mv ${prefix}.gff3 ${prefix}.gff
        fi

        hmmsearch --tblout vgrG.TIGR01646.1.table.out ${hmm1} ${prefix}.faa > /dev/null
        hmmsearch --tblout vgrG.TIGR03361.1.table.out ${hmm2} ${prefix}.faa > /dev/null

        cat vgrG.TIGR*.out | grep -v '#' | sed 's/\\s.*//g' | sort | uniq | while read vgrGlocus; do
            operon=`grep "^\$vgrGlocus	" ${prefix}_operons.tsv | cut -f 2`
            othergenes=`grep "	\$operon\$" ${prefix}_operons.tsv | cut -f 1`
            
            curout=""
            while read curgene; do
                found=`grep "ID=\$curgene" ${prefix}.gff | grep 'CDS' | cut -f 1,4,5,7,9 | sed 's/ID=.*locus_tag=//g;s/;.*product=/\\t/g;s/;.*gene=/\\t/g;s/;.*//g'`
                firstpart=`echo -e "\$found" | cut -f 1-5 -d '	'`
                productpart=`echo -e "\$found" | cut -f 6 -d '	'`
                genepart=`echo -e "\$found" | cut -f 7 -d '	'`
                found="\$firstpart	\$genepart	\$productpart"
                curout="\${curout}\${curgene}	\$found\\n"
            done < <(echo "\$othergenes")
            
            
            echo -e "vgrG-like \$vgrGlocus	\${operon}:"
            
            strand=`echo -e "\$curout" | head -n 1 | cut -f 5`
            if [[ "\$strand" == "+" ]]; then
                echo -e "\$curout" | head -n-1
            else
                echo -e "\$curout" | head -n-1 | tac
            fi
            
            echo -e ""
        done || true

        rm vgrG.TIGR01646.1.table.out 
        rm vgrG.TIGR03361.1.table.out
    } > "${prefix}_T6SS_vgrG_cluster_list.out"
    """

}