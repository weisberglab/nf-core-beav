/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { preCheckInstall              } from '../subworkflows/local/precheckinstall/main.nf'

include { softwareVersionsToYAML        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText        } from '../subworkflows/local/utils_nfcore_beav_pipeline'
include { fileParser                    } from '../subworkflows/local/file_parser'
include { runBakta                      } from '../subworkflows/local/runbakta/main.nf'
include { postBaktaAnnotations          } from '../subworkflows/local/postbaktaannotations/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BEAV {

    take:
    ch_input // channel: input file or files from --input

    main:
    ch_versions = channel.empty()
    
    // Check the models, databases, and other dependencies
    precheck            = preCheckInstall(channel.fromPath("${params.beav_dir}"))
    ch_bakta_db_dl      = precheck.ch_bakta_db_dl
    ch_antismash_db_dl  = precheck.ch_antismash_db_dl
    df_db               = precheck.df_db
    mf_db               = precheck.mf_db
    gap_db              = precheck.gap_db
    t_db                = precheck.t_db

    //
    // Parse input to separate into genbank and fasta files
    //
    ch_parsed       = fileParser(ch_input)
    ch_bakta_input  = ch_parsed.ch_meta_needsbakta
    ch_meta_genbank = ch_parsed.ch_meta_genbank
    ch_meta_fasta   = ch_parsed.ch_meta_fasta
    ch_meta_gff     = ch_parsed.ch_meta_gff
    ch_meta_faa     = ch_parsed.ch_meta_faa
    ch_meta_ffn     = ch_parsed.ch_meta_ffn
    ch_tiger_tuple = channel.empty()
    ch_tiger_tuple = ch_tiger_tuple.mix(ch_meta_fasta, ch_meta_gff, ch_meta_faa).groupTuple().map{ meta, files -> tuple(meta, files[0], files[1], files[2])}


    runBakta(
        ch_bakta_input,
        ch_bakta_db_dl,
        ch_antismash_db_dl,
        df_db,
        mf_db,
        gap_db,
        t_db
    )

    postBaktaAnnotations(
        ch_meta_genbank,
        ch_meta_gff,
        ch_meta_fasta,
        ch_meta_faa,
        ch_meta_ffn,
        ch_antismash_db_dl,
        df_db,
        mf_db,
        gap_db,
        ch_tiger_tuple,
        t_db
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'beav_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
