include { PRECHECK_INSTALL                      } from '../../../modules/local/fileparser/precheckinstall/main.nf'
include { ANTISMASH_ANTISMASHDOWNLOADDATABASES  } from '../../../modules/nf-core/antismash/antismashdownloaddatabases/main.nf'
include { BAKTA_BAKTADBDOWNLOAD                 } from '../../../modules/nf-core/bakta/baktadbdownload/main.nf'
include { PREPARE_GAPMIND                       } from '../../../modules/local/fileparser/getgapminddatabases/main.nf'
include { TIGER_CLONETIGER                      } from '../../../modules/local/tiger/clonetiger/main.nf'
include { DL_TXSS_PAPBLAST                      } from '../../../modules/local/defensefinder/downloaddfdb/main.nf'
include { DL_MACSY_DB                           } from '../../../modules/local/macsyfinder/downloadmfdb/main.nf'

workflow preCheckInstall{
    take:
    ch_main_path

    main:
    done_ch = channel.of("dummy")
    ch_bakta_db_dl = channel.empty()
    ch_antismash_db_dl = channel.empty()
    df_db = channel.empty()
    mf_db = channel.empty()
    gap_db = channel.empty()
    t_db   = channel.empty()

    // Check for existence of a couple databases / files, redownload whole thing if any are missing
    if( !file("${params.beav_dir}/agrobacterium_taxa").exists() || !file("${params.beav_dir}/oncogenic_plasmids").exists() || !file("${params.beav_dir}/tiger.patch2").exists()){
        beav_dl = PRECHECK_INSTALL(ch_main_path)
        beav_dl_ch = beav_dl.done
        done_ch     = done_ch.merge(beav_dl_ch)
    } else {
        beav_dl_ch = channel.of("dummy")
    }

    // TIGER clone
    if ( !file("${params.beav_dir}/TIGER-TIGER2.1/").exists() && !params.skip_tiger ){
        tiger_dl    = TIGER_CLONETIGER(channel.fromPath("${params.beav_dir}/tiger.patch2"), beav_dl_ch)
        t_db = tiger_dl.tiger
    }

    // Defensefinder database
    if ( !file("${params.beav_dir}/models/macsy-finder-models/TXSScan/metadata.yml").exists() && !params.skip_macsyfinder ){
        macsy_dl = DL_MACSY_DB(beav_dl_ch)
        mf_db   = macsy_dl.db
    }

    // Macsyfinder database
    if ( !file("${params.beav_dir}/models/defense-finder-models/defense-finder-models/metadata.yml").exists() && !params.skip_defensefinder )  {
        txss_dl = DL_TXSS_PAPBLAST(beav_dl_ch)
        df_db   = txss_dl.db
    }

    // BAKTA DB
    if( !file("${params.beav_dir}/databases/bakta_db/version.json").exists() && !params.bakta_db){
        println "Downloading ${params.bakta_db_type} bakta database. To change this, run again with --bakta_db_type [light|full]"
        bakta_dl_ch = BAKTA_BAKTADBDOWNLOAD(beav_dl_ch)
        ch_bakta_db_dl = bakta_dl_ch.db
    }

    // ANTISMASH DB
    if( !file("${params.beav_dir}/databases/antismash_db/as-js/0.16/antismash.js").exists() && !params.skip_antismash && !params.antismash_db){
        antismash_dl_ch = ANTISMASH_ANTISMASHDOWNLOADDATABASES(beav_dl_ch)
        ch_antismash_db_dl = antismash_dl_ch.database
    }

    // GAPMIND DB and singularity image creation
    if ( !file("${params.beav_dir}/databases/PaperBLAST_db/README.md").exists() && !params.skip_gapmind ){
        gapprep = PREPARE_GAPMIND(beav_dl_ch)
        gap_db = gapprep.papb
    }

    emit:
    ch_bakta_db_dl
    ch_antismash_db_dl
    df_db
    mf_db
    gap_db
    done_ch
    t_db

}