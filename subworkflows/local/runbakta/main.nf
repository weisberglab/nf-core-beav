#!/usr/bin/env nextflow

include { BAKTA_BAKTA                   } from '../../../modules/nf-core/bakta/bakta/main.nf'

include { postBaktaAnnotations          } from '../postbaktaannotations/main.nf'

workflow runBakta{

    take:
    ch_meta_fasta
    ch_bakta_db_dl
    ch_antismash_db_dl
    df_db
    mf_db
    gap_db
    t_db

    main:    
    // Run BAKTA
    bakta_out = BAKTA_BAKTA(
        ch_meta_fasta,
        params.bakta_db ? ch_meta_fasta.map{ f -> file(params.bakta_db)} : file("${params.beav_dir}/databases/bakta_db").exists() ? ch_meta_fasta.map{ f -> file("${params.beav_dir}/databases/bakta_db")} : ch_bakta_db_dl.combine(ch_meta_fasta).map{ db, meta, f -> file("${params.beav_dir}/databases/bakta_db")},
        [],
        [],
        [],
        []
    )

    // Define output channels
    ch_bakta_meta_gbff = bakta_out.gbff
    ch_bakta_meta_gff  = bakta_out.gff
    ch_bakta_meta_fna  = bakta_out.fna
    ch_bakta_meta_faa  = bakta_out.faa
    ch_bakta_meta_ffn  = bakta_out.ffn
    ch_tiger_tuple     = bakta_out.tiger


    postBaktaAnnotations(
        ch_bakta_meta_gbff,
        ch_bakta_meta_gff,
        ch_bakta_meta_fna,
        ch_bakta_meta_faa,
        ch_bakta_meta_ffn,
        ch_antismash_db_dl,
        df_db,
        mf_db,
        gap_db,
        ch_tiger_tuple,
        t_db
    )


}