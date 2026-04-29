process GAPMIND_GAPMIND{
    tag "${meta.id}"
    label 'process_medium'

    container "oras://community.wave.seqera.io/library/diamond_hmmer_perl-cgi_perl-dbd-sqlite_pruned:35ce1c65a30cd549"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(faa)
    path(db)

    output:
    tuple val(meta), path("${meta.id}_combined_GapMind_results.tab"), emit: gaptab

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    cp -rL ${db} PaperBLAST

    mkdir ${prefix}_gapmind_results
    cd ${prefix}_gapmind_results

    mkdir -p gapmind_out
    mkdir -p GapMind

    echo "Running GapMind amino acid synthesis search for ${prefix}"

    # Run DIAMOND blastp search against amino acid synthesis database
    diamond blastp \\
        --query ../${faa} \\
        --db ../PaperBLAST/tmp/path.aa/curated.faa.dmnd \\
        --out gapmind_out/${prefix}_aa_blast.tsv \\
        --outfmt 6 \\
        --evalue 1e-5 \\
        --max-target-seqs 5 \\
        --threads ${task.cpus}

    echo "Running Gapmind carbon search for ${prefix}"

    # Run DIAMOND blastp search against carbon database
    diamond blastp \\
        --query ../${faa} \\
        --db ../PaperBLAST/tmp/path.carbon/curated.faa.dmnd \\
        --out gapmind_out/${prefix}_carbon_matches.tsv \\
        --evalue 1e-5 --max-target-seqs 1 --outfmt 6

    perl ../PaperBLAST/bin/buildorgs.pl \\
    -out GapMind/orgs \\
    -orgs "file:../${faa}:${prefix}"

    # Build diamond DB for reverse search
    diamond makedb --in GapMind/orgs.faa -d GapMind/orgs.faa.dmnd

    echo "Running forward search(AA)"

    # forward search (amino acid biosynthesis)
    perl ../PaperBLAST/bin/gapsearch.pl \\
        -diamond \\
        -orgs GapMind/orgs \\
        -set aa \\
        -out GapMind/aa.hits \\
        -nCPU ${task.cpus} || true
    
    if [[ ! -f GapMind/aa.hits ]]; then
        echo "Error: aa.hits file not created after gapsearch.pl"
        exit 1
    fi

    echo "Running reverse search(AA)"

    # reverse search
    perl ../PaperBLAST/bin/gaprevsearch.pl \\
        -diamond \\
        -orgs GapMind/orgs \\
        -hits GapMind/aa.hits \\
        -curated ../PaperBLAST/tmp/path.aa/curated.faa.dmnd \\
        -out GapMind/aa.revhits \\
        -nCPU ${task.cpus}

    echo "Making summary(AA)"

    perl ../PaperBLAST/bin/gapsummary.pl \\
        -orgs GapMind/orgs \\
        -set aa \\
        -hits GapMind/aa.hits \\
        -rev GapMind/aa.revhits \\
        -out GapMind/aa.sum

    echo "Checking requirements(AA)"

    # temporarily move the files to here so the perl script can find them

    mv ../PaperBLAST/tmp/path.aa .

    perl ../PaperBLAST/bin/checkGapRequirements.pl \\
        -org GapMind \\
        -results . \\
        -set aa \\
        -out GapMind/aa.sum.warn || true

    mv ./path.aa ../PaperBLAST/tmp

    echo "Making table(AA)"

    awk '\$6 > 0' GapMind/aa.sum.steps | \\
    cut -f 4- | cut -f 1,2,4,5 | grep -v '      ' | \\
    sort -k3 -r | \\
    awk -F"\\t" '!seen[\$1, \$4]++' | \\
    sed 's/\\(\\S\\+\\)\\t\\(\\S\\+\\)\\t\\(\\S\\+\\)\$/\\3:\\2:\\1/g' | \\
    awk 'BEGIN{FS="\\t"; OFS=FS}; { arr[\$2] = arr[\$2] == ""? \$1 : arr[\$2] "," \$1 } END {for (i in arr) print i, arr[i] }' | \\
    sed 's/:/\\t/g' | sort -k1 | grep -v '^locusId' | while read line; do
        locus=`echo -e "\$line" | cut -f 1`
        confidence=`echo -e "\$line" | cut -f 2 | sed 's/0/low/g;s/1/medium/g;s/2/high/g'`
        geneid=`echo -e "\$line" | cut -f 3`
        aa=`echo -e "\$line" | cut -f 4`
        echo -e "\$locus\\tGapMind:\$aa biosynthesis; gene=\$geneid; \$confidence confidence" >> ../${prefix}_combined_GapMind_results.tab
    done

    echo "Running forward search(Carbon)"

    perl ../PaperBLAST/bin/gapsearch.pl \\
        -diamond \\
        -orgs GapMind/orgs \\
        -set carbon \\
        -out GapMind/carbon.hits \\
        -nCPU ${task.cpus} || true
    
    if [[ ! -f GapMind/carbon.hits ]]; then
        echo "Error: carbon.hits file not created after gapsearch.pl"
        exit 1
    fi

    echo "Running reverse search(Carbon)"

    perl ../PaperBLAST/bin/gaprevsearch.pl \\
        -diamond \\
        -orgs GapMind/orgs \\
        -hits GapMind/carbon.hits \\
        -curated ../PaperBLAST/tmp/path.carbon/curated.faa.dmnd \\
        -out GapMind/carbon.revhits \\
        -nCPU ${task.cpus}

    echo "Making summary(Carbon)"

    perl ../PaperBLAST/bin/gapsummary.pl \\
        -orgs GapMind/orgs \\
        -set carbon \\
        -hits GapMind/carbon.hits \\
        -rev GapMind/carbon.revhits \\
        -out GapMind/carbon.sum

    # temporarily move the files to the bin so the perl script can find them

    mv ../PaperBLAST/tmp/path.carbon .

    echo "Checking requirements(carbon)"

    perl ../PaperBLAST/bin/checkGapRequirements.pl \\
        -org GapMind \\
        -results . \\
        -set carbon \\
        -out GapMind/carbon.sum.warn || true

    mv ./path.carbon ../PaperBLAST/tmp

    echo "Making table(Carbon)"

    awk '\$6 > 0' GapMind/carbon.sum.steps | cut -f 4- | cut -f 1,2,4,5 | grep -v '\\t\\t' | sort -k3 -r | \\
    awk -F"\\t" '!seen[\$1, \$4]++' | \\
    sed 's/\\(\\S\\+\\)\\t\\(\\S\\+\\)\\t\\(\\S\\+\\)\$/\\3:\\2:\\1/g' | \\
    awk 'BEGIN{FS="\\t"; OFS=FS}; { arr[\$2] = arr[\$2] == ""? \$1 : arr[\$2] "," \$1 } END {for (i in arr) print i, arr[i] }' | \\
    sed 's/:/\\t/g' | sort -k1 | grep -v '^locusId' | while read line; do
        locus=`echo -e "\$line" | cut -f 1`
        confidence=`echo -e "\$line" | cut -f 2 | sed 's/0/low/g;s/1/medium/g;s/2/high/g'`
        geneid=`echo -e "\$line" | cut -f 3`
        carbon=`echo -e "\$line" | cut -f 4`
        echo -e "\$locus\\tGapMind:\$carbon catabolism; gene=\$geneid; \$confidence confidence" >> ../${prefix}_combined_GapMind_results.tab
    done
    cd ..

    echo "All GapMind searches completed for ${prefix}"

    # Cleanup
    rm -rf ./PaperBLAST
    """

}