#!/usr/bin/env nextflow

include { INTEGRONFINDER                    } from '../../../modules/nf-core/integronfinder/main.nf'
include { ANTISMASH_ANTISMASH               } from '../../../modules/nf-core/antismash/antismash/main.nf'
include { ANTISMASH_GENBANKS_TO_TABLE       } from '../../../modules/local/antismash/genbankstotable/main.nf'
include { FASTANI                           } from '../../../modules/nf-core/fastani/main.nf'
include { BLAST_BLASTN                      } from '../../../modules/nf-core/blast/blastn/main.nf'
include { PHISPY                            } from '../../../modules/nf-core/phispy/main.nf'
include { HMMER_HMMSEARCH                   } from '../../../modules/nf-core/hmmer/hmmsearch/'
include { SOURMASH_SKETCH                   } from '../../../modules/nf-core/sourmash/sketch/'
include { SOURMASH_COMPARE                  } from '../../../modules/nf-core/sourmash/compare/main.nf'
include { TIGER_TIGER                       } from '../../../modules/local/tiger/tiger/main.nf'
include { MACSYFINDER_MACSYFINDER           } from '../../../modules/local/macsyfinder/macsyfinder/main.nf'
include { DEFENSEFINDER_DEFENSEFINDER       } from '../../../modules/local/defensefinder/defensefinder/main.nf'
include { OPERONMAPPER_SUBMITOPERONMAPPERJOB} from '../../../modules/local/operonmapper/submitoperonmapperjob/main.nf'
include { FUZZNUC_MAIN_ELEMENTS             } from '../../../modules/local/fuzznuc/fuzznucmainelements/main.nf'
include { GAPMIND_GAPMIND                   } from '../../../modules/local/gapmind/gapmind/main.nf'
include { IDENTIFY_TI_RI_PLASMIDS           } from '../../../modules/local/tiriplasmids/identifytiriplasmids/main.nf'
include { HMMER_NHMMER                      } from '../../../modules/local/hmmer/nhammer/main.nf'
include { COMBINE_ANNOTATIONS               } from '../../../modules/local/combineannotations/combineannotations/main.nf'
include { CIRCOS_CIRCOS                     } from '../../../modules/local/circos/circos/main.nf'
include { IDENTIFY_VGRG_CLUSTERS            } from '../../../modules/local/vgrgclusters/identify/main.nf'
include { MOBSUITE_RECON                    } from '../../../modules/nf-core/mobsuite/recon/main.nf'

workflow postBaktaAnnotations{

    take:
    ch_meta_gbff
    ch_meta_gff
    ch_meta_fna
    ch_meta_faa
    ch_meta_ffn
    ch_antismash_db_dl
    df_db
    mf_db
    gap_db
    ch_tiger_tuple
    t_db

    main:
    combined_ch = ch_meta_gbff
    combined_ch = combined_ch.mix(ch_meta_faa)
    combined_ch = combined_ch.mix(ch_meta_ffn)

    //
    // NF-CORE MODULE: MOBSUITE_RECON || Plasmid Characterization
    //
    if(!params.skip_mobsuite){
        mobsuiteresults = MOBSUITE_RECON(
            ch_meta_fna
        )
        combined_ch = combined_ch.mix(mobsuiteresults.contig_report)
        combined_ch = combined_ch.mix(mobsuiteresults.biomarkers)
    }

    //
    // LOCAL MODULE: MACSYFINDER    || Secretion Systems
    //
    if(!params.skip_macsyfinder){

        macsyresults = MACSYFINDER_MACSYFINDER(
            ch_meta_faa,
            file("${params.beav_dir}/models/macsy-finder-models/TXSScan/metadata.yml").exists() ? Channel.fromPath("${params.beav_dir}/models/macsy-finder-models").combine(ch_meta_faa).map{ d, m, f -> d } : mf_db.combine(ch_meta_faa).map{ d, m, f -> file("${params.beav_dir}/models/macsy-finder-models") }
        )

        combined_ch = combined_ch.mix(macsyresults.tsv)

    }

    //
    // LOCAL MODULE: DEFENSEFINDER  || Defense Systems
    //
    if(!params.skip_defensefinder){
        defenseout = DEFENSEFINDER_DEFENSEFINDER(
            ch_meta_faa,
            file("${params.beav_dir}/models/defense-finder-models/defense-finder-models/metadata.yml").exists() ? Channel.fromPath("${params.beav_dir}/models/defense-finder-models").combine(ch_meta_faa).map{ d, m, f -> d } : df_db.combine(ch_meta_faa).map{ d, m, f -> file("${params.beav_dir}/models/defense-finder-models") }
        )
        
        combined_ch = combined_ch.mix(defenseout.tsv)
        
    }

    //
    // NF-CORE MODULE: ANTISMASH    || Secondary Metabolite Cluster Prediction
    //
    if (!params.skip_antismash){

        antismash_out = ANTISMASH_ANTISMASH(
            ch_meta_gbff,
            params.antismash_db ? ch_meta_gbff.map{f -> file(params.antismash_db)} : file("${params.beav_dir}/databases/antismash_db").exists() ? ch_meta_gbff.map{f -> file("${params.beav_dir}/databases/antismash_db")} : ch_antismash_db_dl.combine(ch_meta_gbff).map{ db, m, file -> db},
            []
        )

        antismash_final = ANTISMASH_GENBANKS_TO_TABLE(
            antismash_out.gbkresultsdir
        )

        combined_ch = combined_ch.mix(antismash_final.subset)
    }

    //
    // LOCAL MODULES: OPERONMAPPER, IDENTIFY_VGRG_CLUSTERS       || Operon Prediction and vgrG Clusters
    //
    if (params.operon_email){
        ch_email = Channel.of("${params.operon_email}").combine(ch_meta_fna).map{ e, m, f -> e }

        // Join gff to fna with basename
        ch_fna_gff = Channel.empty()
        ch_fna_gff = ch_fna_gff.mix(ch_meta_fna, ch_meta_gff).groupTuple()

        operonmapperout = OPERONMAPPER_SUBMITOPERONMAPPERJOB(
            ch_fna_gff,
            ch_email
        )

        ch_operons = operonmapperout.tsv
        ch_vgrg_input = ch_operons.mix(ch_meta_faa, ch_meta_gff).groupTuple()

        vgrgout = IDENTIFY_VGRG_CLUSTERS(
            ch_vgrg_input,
            ch_vgrg_input.map{f -> file("${params.beav_dir}/models/vgrgmodels/TIGR01646.1.hmm")},
            ch_vgrg_input.map{f -> file("${params.beav_dir}/models/vgrgmodels/TIGR03361.1.hmm")}
            )

        combined_ch = combined_ch.mix(vgrgout.clist)

    } else {
        println "Operon email not provided. Skipping operon mapper."
    }

    //
    // LOCAL MODULE: GAPMIND            || AA Biosynthesis Pathways
    //
    if (!params.skip_gapmind){

        gapmindout = GAPMIND_GAPMIND(
            ch_meta_faa,
            file("${params.beav_dir}/databases/PaperBLAST_db/tmp/path.aa/curated.faa").exists() ? ch_meta_faa.map{f -> file("${params.beav_dir}/databases/PaperBLAST_db")} : gap_db.combine(ch_meta_faa).map{ env, m, file -> env},
        )

        combined_ch = combined_ch.mix(gapmindout.gaptab)

    }

    //
    // LOCAL MODULE: NHMMER             || Dif sites using Hidden Markov Models
    //
    if (!params.skip_hmmer){
        // Popluate input HMM channels
        ch_tDNA_leftborder  = []
        ch_tDNA_rightborder = []
        ch_overdrive        = []
        if (params.agrobacterium || params.agro){
            ch_tDNA_leftborder  = Channel.fromPath("${params.beav_dir}/models/agromodels/T-DNA_leftborder.hmm").combine(ch_meta_fna).map{ t, m, f -> t }
            ch_tDNA_rightborder = Channel.fromPath("${params.beav_dir}/models/agromodels/T-DNA_rightborder.hmm").combine(ch_meta_fna).map{ t, m, f -> t }
            ch_overdrive        = Channel.fromPath("${params.beav_dir}/models/agromodels/overdrive.hmm").combine(ch_meta_fna).map{ t, m, f -> t }
        }

        ch_nodbox = Channel.fromPath("${params.beav_dir}/models/agromodels/nodbox.hmm").combine(ch_meta_fna).map{ t, m, f -> t }
        ch_ttsbox = Channel.fromPath("${params.beav_dir}/models/agromodels/ttsbox.hmm").combine(ch_meta_fna).map{ t, m, f -> t }


        hmmerout = HMMER_NHMMER(
            ch_meta_fna,
            ch_nodbox,
            ch_ttsbox,
            ch_tDNA_leftborder,
            ch_tDNA_rightborder,
            ch_overdrive
        )

        combined_ch = combined_ch.mix(hmmerout.tsv)

    }

    //
    // LOCAL MODULE: FUZZNUC            || pip, tts, nod, hrp boxes
    // 
    if (!params.skip_fuzznuc){

        if (params.agrobacterium){
            ch_agro_true = Channel.of(true).combine(ch_meta_fna).map{ t, m, f -> t }

            fuzzoutm = FUZZNUC_MAIN_ELEMENTS(
                ch_meta_fna,
                ch_agro_true
            )
            
            combined_ch = combined_ch.mix(fuzzoutm.tsv)

        } else {
            
            ch_agro_false = Channel.of(false).combine(ch_meta_fna).map{ f, m, fi -> f }

            fuzzoutm = FUZZNUC_MAIN_ELEMENTS(
                ch_meta_fna,
                ch_agro_false
            )
            
            combined_ch = combined_ch.mix(fuzzoutm.tsv)
        }
        
    }

    //
    // NF-CORE MODULE: BLASTN           || Nucleotide Sequence Similarity Search
    //
    if (!params.skip_blastn){
        // Populate input db channels
        ch_oriT = Channel.fromPath("${params.beav_dir}/oriT_db.fna").combine(ch_meta_fna).map{ t, m, f -> t }
        ch_dif  = Channel.fromPath("${params.beav_dir}/dif_db.fna").combine(ch_meta_fna).map{ t, m, f -> t }
        ch_blastn_false = Channel.of(false).combine(ch_meta_fna).map{ f, m, fi -> f }

        blout = BLAST_BLASTN(
            ch_meta_fna,
            ch_oriT,
            ch_dif,
            [],
            [],
            ch_blastn_false
        )

        combined_ch = combined_ch.mix(blout.oriTtable)
        combined_ch = combined_ch.mix(blout.diftable)
    }

    //
    // NF-CORE MODULES: SOURMASH_SKETCH, SOURMASH_COMPARE   || Genome Sketching and Comparison
    //
    if(!params.skip_sourmash){
        sourmashout = SOURMASH_SKETCH(
                                ch_meta_fna
                             )

        signatures_ch = sourmashout.signatures

        ch_true = Channel.of(true).combine(ch_meta_fna).map{ t, m, f -> t }

        smashfout = SOURMASH_COMPARE(
            signatures_ch,
            [],
            ch_true,
            ch_true
        )

        combined_ch = combined_ch.mix(smashfout.csv)
    }

    
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MOBILE GENETIC ELEMENTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // 
    // LOCAL MODULE: TIGER              || ICEs
    //
    if(params.tiger_blast_db && !params.skip_tiger){
        
        ch_tiger_dir = file("${params.beav_dir}/TIGER-TIGER2.1").exists() ? Channel.fromPath("${params.beav_dir}/TIGER-TIGER2.1").combine(ch_meta_fna).map{ db, _m, _file -> db } : t_db.combine(ch_meta_fna).map { db, _m, _f -> file("${params.beav_dir}/TIGER-TIGER2.1")}
        ch_tiger_blast_db = channel.fromPath("${params.tiger_blast_db}").combine(ch_meta_fna).map{ db, _m, _file -> db }

        tigerout = TIGER_TIGER( 
            ch_tiger_tuple,
            ch_tiger_dir,
            ch_tiger_blast_db
        )

        combined_ch = combined_ch.mix(tigerout.table)
                
    } else if(!params.skip_tiger) {
        println("TIGER BLAST database not specified. A database is required for TIGER ICE analysis.")
        println("Specify with --tiger_blast_db path/to/database")
        println("If you do not want to run TIGER you may safely ignore this message.")
        println("To prevent future messages from appearing, specify --skip_tiger")
    }


    //
    // NF-CORE MODULE: INTEGRONFINDER   || Integron Finder
    //
    if(!params.skip_integronfinder){
        integronout = INTEGRONFINDER( ch_meta_fna )

        combined_ch = combined_ch.mix(integronout.summary)
        combined_ch = combined_ch.mix(integronout.integrons)

    }

    //
    // NF-CORE MODULE: PHISPY           || Prophage Identification
    //
    if(!params.skip_phispy){
        phispyout = PHISPY( ch_meta_gbff )

        combined_ch = combined_ch.mix(phispyout.coordinates)

    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        AGROBACTERIUM SPECIFIC PIPELINE
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if(params.agrobacterium || params.agro){

        //
        // NF-CORE MODULE: fastANI      || Whole-genome Average Nucleotide Identity
        //
        ch_refs = Channel.fromPath("${params.beav_dir}/agrobacterium_taxa/*.fna").combine(ch_meta_fna).map{ r, m, f -> r }
        fastaniout = FASTANI(
            ch_meta_fna,    // query genomes
            ch_refs  // reference genomes
        )

        combined_ch = combined_ch.mix(fastaniout.ani)

        //
        // LOCAL MODULE: IDENTIFY TI RI PLASMIDS
        //
        ch_oncgs = Channel.fromPath("${params.beav_dir}/oncogenic_plasmids").combine(ch_meta_fna).map{ o, m, f -> o }
        tiriout = IDENTIFY_TI_RI_PLASMIDS(
            ch_meta_fna,
            ch_oncgs
        )

        contigs_ch  = Channel.empty()
        contigs_ch  = tiriout.contiglist
        
    } else {
        contigs_ch  = Channel.empty()
    }

    combined_annotations = COMBINE_ANNOTATIONS(
        combined_ch.groupTuple() // Group by meta id
    )

    circos_ch = Channel.empty()
    circos_ch = circos_ch.mix(combined_annotations.finalgbk, combined_annotations.finalffn, combined_annotations.finalfaa, contigs_ch).groupTuple()

    CIRCOS_CIRCOS(
        circos_ch
    )

}
