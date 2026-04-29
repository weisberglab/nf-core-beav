process PREPARE_GAPMIND{
    label 'process_medium'

    publishDir "${params.beav_dir}/databases", mode: 'copy'

    container "oras://community.wave.seqera.io/library/diamond_perl-dbd-sqlite_curl_perl-dbi_pruned:597d98710e7c495e"
    conda "${moduleDir}/environment.yml"

    input:
    val(dummy)

    output:
    path("PaperBLAST_db"),  emit: papb

    script:
    """
    # Make amino acid path

    curl -L https://github.com/morgannprice/PaperBLAST/archive/refs/heads/master.zip -o PaperBLAST.zip
    unzip PaperBLAST.zip
    mv PaperBLAST-master PaperBLAST

    mkdir -p PaperBLAST/tmp/path.aa
    cd PaperBLAST/tmp/path.aa
    curl -O https://papers.genomics.lbl.gov/tmp/path.aa/curated.faa
    curl -O https://papers.genomics.lbl.gov/tmp/path.aa/curated.db
    curl -O https://papers.genomics.lbl.gov/tmp/path.aa/steps.db
    perl ../../bin/extractHmms.pl steps.db .
    
    # Make carbon path
    mkdir -p ../path.carbon
    cd ../path.carbon
    curl -O https://papers.genomics.lbl.gov/tmp/path.carbon/curated.faa
    curl -O https://papers.genomics.lbl.gov/tmp/path.carbon/curated.db
    curl -O https://papers.genomics.lbl.gov/tmp/path.carbon/steps.db
    perl ../../bin/extractHmms.pl steps.db .
    cd ../../..

    # Other dependencies
    mkdir -p PaperBLAST/fbrowse_data
    mkdir -p PaperBLAST/private
    mkdir -p PaperBLAST/tmp/downloaded

    # Make diamond
    # wget -q https://github.com/bbuchfink/diamond/releases/download/v2.0.15/diamond-linux64.tar.gz
    # tar xzf diamond-linux64.tar.gz
    # rm diamond-linux64.tar.gz

    # Prebuild diamond databases(amino acid and carbon)
    diamond makedb --in PaperBLAST/tmp/path.aa/curated.faa \\
               -d PaperBLAST/tmp/path.aa/curated.faa.dmnd
    diamond makedb --in PaperBLAST/tmp/path.carbon/curated.faa \\
               -d PaperBLAST/tmp/path.carbon/curated.faa.dmnd

    # Patch gapsearch.pl to find hmmsearch and diamond in the environment
    sed -i 's|my \$hmmsearch = "\$binDir/hmmsearch";|my \$hmmsearch = `which hmmsearch`; chomp(\$hmmsearch);|' PaperBLAST/bin/gapsearch.pl
    sed -i 's|my \$diamond = "\$binDir/diamond";|my \$diamond = `which diamond`; chomp(\$diamond);|' PaperBLAST/bin/gapsearch.pl
    sed -i 's|my \$diamond = "\$Bin/diamond";|my \$diamond = `which diamond`; chomp(\$diamond);|' PaperBLAST/bin/gaprevsearch.pl

    # Rename it so it can be copied easily to each work dir
    mv PaperBLAST PaperBLAST_db
    """
}