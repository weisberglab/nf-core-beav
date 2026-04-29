process TIGER_CLONETIGER{
    
    publishDir "${params.beav_dir}", mode: 'copy'

    container "oras://community.wave.seqera.io/library/curl_gcc_patch:84930b7881d83ff8"
    conda "${moduleDir}/environment.yml"

    input:
    path(patch)
    val(dummy)

    output:
    path("TIGER-TIGER2.1"), emit: tiger
    val(true), emit: done

    script:
    """
    # Clone and patch tiger
    curl -v -L -O https://github.com/sandialabs/TIGER/archive/refs/tags/TIGER2.1.tar.gz
    tar xzf TIGER2.1.tar.gz --exclude 'TIGER-TIGER2.1/db/Pfam-A.hmm' --exclude 'TIGER-TIGER2.1/bin/aragorn1.2*' --exclude 'TIGER-TIGER2.1/bin/aragorn1.1' --exclude 'TIGER-TIGER2.1/bin/hmmsearch' --exclude 'TIGER-TIGER2.1/bin/pfscan'
    gcc -O3 -ffast-math -finline-functions -o TIGER-TIGER2.1/bin/aragorn1.1 TIGER-TIGER2.1/bin/aragorn1.1.6.c
    rm -f TIGER2.1.tar.gz

    curl -L "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam25.0/Pfam-A.hmm.gz" -o "TIGER-TIGER2.1/db/Pfam-A.hmm.gz"
    gunzip TIGER-TIGER2.1/db/Pfam-A.hmm.gz
    rm -f TIGER-TIGER2.1/db/Pfam-A.hmm.gz

    patch -p0 --batch < ${patch}

    find TIGER-TIGER2.1/bin -name "*.pl" -type f -exec sed -i.bak -E 's/&>\\s*([^\\s]+)/> \\1/g' {} +    

    chmod +x TIGER-TIGER2.1/bin/*pl
    """
}