process BEDTOOLS {
    tag "$meta.id"
    label 'process_low'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2c/2cb27be9b55640ce5717c51a6a6bf79b7945d06eab425f0616d46251b4236e9e/data' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(oriT), path(dif)

    output:
    tuple val(meta), path("${prefix}_oriT.table"), emit: oriT
    tuple val(meta), path("${prefix}_dif.table") , emit: dif
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version"), topic: versions, emit: versions_bedtools

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # OriT regions
    touch ${prefix}_oriT.table
    touch oriT_hits.stranded

    cat ${oriT} | awk '\$4 >= 20' | awk '\$11 <= 0.1' | while read line; do
        contig=`echo -e "\$line" | cut -f 2`
        queryname=`echo -e "\$line" | cut -f 1`
        startpos=`echo -e "\$line" | cut -f 9`
        endpos=`echo -e "\$line" | cut -f 10`
        evalue=`echo -e "\$line" | cut -f 11`
        if [[ "\$endpos" -gt "\$startpos" ]]; then
            strand="+1"
        else
            strand="-1"
            temppos="\$startpos"
            startpos="\$endpos"
            endpos="\$temppos"
        fi
        echo -e "\$contig	\$startpos	\$endpos	\$strand	\$queryname	\$evalue" >> oriT_hits.stranded 
    done 

    if [ -s oriT_hits.stranded ]; then
        cat oriT_hits.stranded | \\
            sort -k1,1 -k2,2g | \\
            bedtools merge -c 4,5,6 -o collapse,collapse,collapse | \\
            while read curborderline; do
                contig=`echo -e "\$curborderline" | cut -f 1`
                startpos=`echo -e "\$curborderline" | cut -f 2`
                endpos=`echo -e "\$curborderline" | cut -f 3`
                
                oriTevalues=`echo -e "\$curborderline" | cut -f 6`
                oriTnames=`echo -e "\$curborderline" | cut -f 5`
                oriTstrand=`echo -e "\$curborderline" | cut -f 4`
                
                evallist=`echo -e "\$oriTevalues" | sed 's/,/\\n/g'`
                namelist=`echo -e "\$oriTnames" | sed 's/,/\\n/g'`
                strandlist=`echo -e "\$oriTstrand" |  tr ',' '\\n'`
                
                newname=`echo "\$namelist" | paste - <(echo "\$evallist") | paste - <(echo "\$strandlist") | sort -k2,2g | head -n1 `
                newstrand=`echo -e "\$newname" | cut -f 3`
                newannot=`echo -e "\$newname" | cut -f 1-2 | tr '\t' ';'`
                echo -e "\$contig	oriT	\$startpos	\$endpos	\$newstrand	\$newannot" >> ${prefix}_oriT.table
            done
    fi

    # Dif Sites
    touch ${prefix}_dif.table
    touch dif_hits.stranded

    cat ${prefix}_difSites.txt | awk '\$4 >= 20' | awk '\$11 <= 0.01' | while read line; do
        contig=`echo -e "\$line" | cut -f 2`
        queryname=`echo -e "\$line" | cut -f 1`
        startpos=`echo -e "\$line" | cut -f 9`
        endpos=`echo -e "\$line" | cut -f 10`
        evalue=`echo -e "\$line" | cut -f 11`
        if [[ "\$endpos" -gt "\$startpos" ]]; then
            strand="+1"
        else
            strand="-1"
            temppos="\$startpos"
            startpos="\$endpos"
            endpos="\$temppos"
        fi
        echo -e "\$contig	\$startpos	\$endpos	\$strand	\$queryname	\$evalue" >> dif_hits.stranded 
    done             

    if [ -s dif_hits.stranded ]; then
        cat dif_hits.stranded | \\
            sort -k1,1 -k2,2g | \\
            bedtools merge -c 4,5,6 -o collapse,collapse,collapse | \\
            while read curborderline; do
                contig=`echo -e "\$curborderline" | cut -f 1`
                startpos=`echo -e "\$curborderline" | cut -f 2`
                endpos=`echo -e "\$curborderline" | cut -f 3`
                
                difevalues=`echo -e "\$curborderline" | cut -f 6`
                difnames=`echo -e "\$curborderline" | cut -f 5`
                difstrand=`echo -e "\$curborderline" | cut -f 4`
                
                evallist=`echo -e "\$difevalues" | sed 's/,/\\n/g'`
                namelist=`echo -e "\$difnames" | sed 's/,/\\n/g'`
                strandlist=`echo -e "\$difstrand" |  tr ',' '\\n'`
                
                newname=`echo "\$namelist" | paste - <(echo "\$evallist") | paste - <(echo "\$strandlist") | sort -k2,2g | head -n1 `
                newstrand=`echo -e "\$newname" | cut -f 3`
                newannot=`echo -e "\$newname" | cut -f 1-2 | tr '\t' ';'`
                echo -e "\$contig	\$startpos	\$endpos	100	dif	\$newstrand">> ${prefix}_dif.table
            done
    fi

    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo $args
    
    touch ${prefix}_oriT.table
    touch ${prefix}_dif.table
    """
}
