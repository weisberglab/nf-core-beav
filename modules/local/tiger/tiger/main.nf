process TIGER_TIGER {
    tag "${meta.id}"
    label 'process_medium'

    container 'oras://community.wave.seqera.io/library/aragorn_bedtools_blast_hmmer_pruned:419f84dd344e8c0c'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fna), path(gff3), path(faa)
    path(tigerdb)
    path(blastdb)

    output:
    tuple val(meta), path("${prefix}_TIGER2_final.table.out") , emit: table

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    # Prepare tiger inputs n stuff

    if [[ -f "${prefix}.gff" ]]; then
        cp ${prefix}.gff ${prefix}.gff3
    fi

    mkdir -p TIGER2/protein
    cp ${prefix}.fna TIGER2/genome.fa
    cp ${prefix}.faa TIGER2/protein/protein.faa
    # Remove > and < from the gff
    awk 'BEGIN{FS=OFS="\\t"} { gsub(/[<>]/, "", \$4); gsub(/[<>]/, "", \$5); print }' ${prefix}.gff3 > TIGER2/protein/protein.gff
    cd TIGER2

    # Run islander
    perl ../${tigerdb}/bin/islander.pl -verbose -tax B -cpu ${task.cpus} genome.fa

    # Run tiger
    perl ../${tigerdb}/bin/tiger.pl -verbose \\
        -db ${blastdb} \\
        -fasta genome.fa \\
        -cpu ${task.cpus} \\
        -outDir ${prefix}

    # Run typing
    perl ../${tigerdb}/bin/typing.pl genome.island.nonoverlap.gff

    # Run resolve
    perl ../${tigerdb}/bin/resolve.pl mixed lenient

    # Run final typing
    perl ../${tigerdb}/bin/typing.pl resolve3.gff

    echo "Onto parsing output."

    # Parse output
    set +u
    declare -A aatable
    aatable=( ["C"]="Cys" 
            ["D"]="Asp"
            ["S"]="Ser"
            ["Q"]="Gln"
            ["K"]="Lys"
            ["I"]="Ile"
            ["P"]="Pro"
            ["T"]="Thr"
            ["F"]="Phe"
            ["N"]="Asn"
            ["G"]="Gly"
            ["H"]="His"
            ["L"]="Leu"
            ["R"]="Arg"
            ["W"]="Trp"
            ["A"]="Ala"
            ["V"]="Val"
            ["E"]="Glu"
            ["Y"]="Tyr"
            ["M"]="Met" )
            
    cd ..

    # prepare prefix.nofasta.gff3
    awk -F'\t' 'BEGIN{OFS="\t"}
    /^##FASTA/ {exit}           # stop at FASTA section
    !/region\t[0-9]/ && /CDS/ && !/remark/ {print}' \
    ${prefix}.gff3 \
    > ${prefix}.nofasta_pre.gff3
    awk 'BEGIN{FS=OFS="\\t"} { gsub(/[<>]/, "", \$4); gsub(/[<>]/, "", \$5); print }' ${prefix}.nofasta_pre.gff3 > ${prefix}.nofasta.gff3

    touch ${prefix}_TIGER2_final.table.out 
    while read line; do
        replicon=`echo -e "\$line" |  cut -f 1`
        islandtype=`echo -e "\$line" |  cut -f 2`
        startpos=`echo -e "\$line" |  cut -f 4`
        endpos=`echo -e "\$line" |  cut -f 5`
        
        tigerID=`echo -e "\$line" |  cut -f 9 | sed 's/^.*ID=//g;s/;.*//g'`

        target=`echo -e "\$line" |  cut -f 9 | sed 's/;coord=.*//g;s/^.*target=//g'`
        targetlen=`expr length "\$target"`
        if [[ "\$targetlen" == 1 ]]; then
            target="tRNA_\${aatable[\$target]}"
        fi
            if [ \$islandtype == "TIGER" ] || [ \$islandtype == "Islander,TIGER" ]; then
                leftborderseq=`echo -e "\$line" |  cut -f 9 | sed 's/^.*isleLseq=//g;s/;.*//g' | sed 's/[A-Z]//g'`
                rightborderseq=`echo -e "\$line" |  cut -f 9 | sed 's/^.*isleRseq=//g;s/;.*//g' | sed 's/[A-Z]//g'`
                targetseq=`echo -e "\$line" |  cut -f 9 | sed 's/^.*unintSeq=//g;s/;.*//g' | sed 's/[A-Z]//g'`
            fi
        icelength=`echo -e "\$line" |  cut -f 9 | sed 's/^.*;len=//g;s/;.*//g'`
        tandem=`echo -e "\$line" |  cut -f 9 | sed 's/^.*;tandem=//g;s/;.*//g' | sed 's/1,1/no/g;s/\\([0-9]\\),/count=\\1,/g;s/,/,pos=/g'`
        
        #get whether ICE or Phage or PhageICE
        MGEtype=`echo -e "\$line" | cut -f 1-5 | grep -f - ./TIGER2/islesFinal.gff |  cut -f 9 | sed 's/^.*;type=//g;s/;.*//g;s/[0-9]\\+\$//g' | sed 's/other/IME/g'`
        
        #if phage, get phage loci
        if [[ \$MGEtype == *"Phage"*  ]]; then
            phagecomp=`cat ./TIGER2/Isles/\$tigerID/phage.txt | grep -v '#nick' | cut -f 5- | sed 's/\\t/,/g' | sed 's/^,\\+//g;s/,\\+\$//g'`
            phageloci="phage_loci=\$phagecomp;"	
        else
            phageloci=""
        fi

        #get loci for all integrases
        integraseposlist=`echo -e "\$line" |  cut -f 9 | sed 's/^.*;ints=//g;s/;.*//g'`
        curloci=""
        while read intlocus; do
            inttype=`echo -e "\$intlocus" | sed 's/\\..*//g'` 
            intposstart=`echo -e "\$intlocus" | sed 's/^.*://g;s/-.*//g'` 
            intposend=`echo -e "\$intlocus" | sed 's/.*-//g'` 
            intlocus=`echo "\$replicon	\$intposstart	\$intposend" | bedtools intersect -nonamecheck -a stdin -b ./${prefix}.nofasta.gff3 -wb -f 0.90  | grep 'CDS' | grep -v 'remark' | sed 's/^.*locus_tag=//g;s/;.*//g' | tr '\\n' ',' | sed 's/,\$//g'`
            curloci="\$curloci,\$inttype:\$intlocus"
        done < <( echo -e "\$integraseposlist" | tr ',' '\\n' | sort -k2,2g -t : )  
        curloci=`echo -e "\$curloci" | sed 's/^,//g'`
        if [ \$islandtype == "TIGER" ] || [ \$islandtype == "Islander,TIGER" ]; then
                echo -e "\$replicon      \$startpos       \$endpos \$MGEtype;target=\$target;attL=\$leftborderseq;attR=\$rightborderseq;attB=\$targetseq;length=\$icelength;integrase=\$curloci;tandem=\$tandem;\${phageloci}pred_model=\$islandtype" >> ${prefix}_TIGER2_final.table.out
            else
                    echo -e "\$replicon      \$startpos       \$endpos \$MGEtype;target=\$target;length=\$icelength;integrase=\$curloci;tandem=\$tandem;\${phageloci}pred_model=\$islandtype" >> ${prefix}_TIGER2_final.table.out
            fi
    done < ./TIGER2/resolve3.gff
    """
}
