process TIGER_CLONETIGER{
    
    publishDir "${params.beav_dir}", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/42/428eaf918f63e696914ec3236735092b18b0f0c128da16a5cc7ae67dd0b1f44c/data'
        : 'community.wave.seqera.io/library/curl_gcc_patch:3dd8a0cb75380c43'}"

    input:
    path(patch)

    output:
    path("TIGER-TIGER2.1"), emit: tiger

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