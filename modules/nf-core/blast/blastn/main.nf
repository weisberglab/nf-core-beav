process BLAST_BLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/bedtools_blast:e4f7c8c6d31f67f9"

    input:
    tuple val(meta) , path(fasta)
    path oriT
    path dif
    path taxidlist
    val taxids
    val negative_tax

    output:
    tuple val(meta), path("${prefix}_oriT.table")  , emit: oriTtable
    tuple val(meta), path("${prefix}_dif.table")   , emit: diftable
    path "versions.yml"           , emit: versions
    val(true), emit: done

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    def negative_tax_cmd = negative_tax ? "negative_" : ""
    def taxidlist_cmd = taxidlist ? "-${negative_tax_cmd}taxidlist ${taxidlist}" : ""
    def taxids_cmd = taxids ? "-${negative_tax_cmd}taxids ${taxids}" : ""
    if (taxidlist_cmd.any() && taxids_cmd.any()) {
        log.error("ERROR: taxidlist and taxids can not be used at the same time, choose only one argument to use for tax id filtering.")
    }

    """
    makeblastdb -in ${fasta} -dbtype nucl 1> /dev/null

    # OriT regions
    touch ${prefix}_oriT.table
    touch oriT_hits.stranded

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
        
    cat ${prefix}_oriT.txt | awk '\$4 >= 20' | awk '\$11 <= 0.1' | while read line; do
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
    touch ${meta.id}_dif.table
    touch dif_hits.stranded

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}