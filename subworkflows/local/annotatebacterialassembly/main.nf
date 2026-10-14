#!/usr/bin/env nextflow
include { BAKTA_BAKTA                        } from '../../../modules/nf-core/bakta/bakta/main.nf'

include { INTEGRONFINDER                     } from '../../../modules/nf-core/integronfinder/main.nf'
include { ANTISMASH_ANTISMASH                } from '../../../modules/nf-core/antismash/antismash/main.nf'
include { ANTISMASH_GENBANKS_TO_TABLE        } from '../../../modules/local/antismash/genbankstotable/main.nf'
include { FASTANI                            } from '../../../modules/nf-core/fastani/main.nf'
include { BLAST_BLASTN                       } from '../../../modules/nf-core/blast/blastn/main.nf'
include { PHISPY                             } from '../../../modules/nf-core/phispy/main.nf'
include { BEDTOOLS                           } from '../../../modules/local/bedtools/main.nf'
include { HMMER_HMMSEARCH                    } from '../../../modules/nf-core/hmmer/hmmsearch/'
include { SOURMASH_SKETCH                    } from '../../../modules/nf-core/sourmash/sketch/'
include { SOURMASH_COMPARE                   } from '../../../modules/nf-core/sourmash/compare/main.nf'
include { TIGER_TIGER                        } from '../../../modules/local/tiger/tiger/main.nf'
include { MACSYFINDER_MACSYFINDER            } from '../../../modules/local/macsyfinder/macsyfinder/main.nf'
include { DEFENSEFINDER_DEFENSEFINDER        } from '../../../modules/local/defensefinder/defensefinder/main.nf'
include { OPERONMAPPER_SUBMITOPERONMAPPERJOB } from '../../../modules/local/operonmapper/submitoperonmapperjob/main.nf'
include { FUZZNUC_MAIN_ELEMENTS              } from '../../../modules/local/fuzznuc/fuzznucmainelements/main.nf'
include { GAPMIND_GAPMIND                    } from '../../../modules/local/gapmind/gapmind/main.nf'
include { IDENTIFY_TI_RI_PLASMIDS            } from '../../../modules/local/tiriplasmids/identifytiriplasmids/main.nf'
include { HMMER_NHMMER                       } from '../../../modules/local/hmmer/nhammer/main.nf'
include { COMBINE_ANNOTATIONS                } from '../../../modules/local/combineannotations/combineannotations/main.nf'
include { CIRCOS_CIRCOS                      } from '../../../modules/local/circos/circos/main.nf'
include { IDENTIFY_VGRG_CLUSTERS             } from '../../../modules/local/vgrgclusters/identify/main.nf'
include { MOBSUITE_RECON                     } from '../../../modules/nf-core/mobsuite/recon/main.nf'

workflow ANNOTATEBACTERIALASSEMBLY{

    take:
    ch_bakta_input
    ch_bakta_db_dl
    ch_meta_gbff
    ch_meta_gff
    ch_meta_fna
    ch_meta_faa
    ch_meta_ffn         // never explicitly used, appears in output files for completeness
    ch_antismash_db_dl
    df_db
    mf_db
    gap_db
    ch_tiger_tuple
    t_db
    ch_filecount

    main:
    combined_ch = ch_meta_gbff
    combined_ch = combined_ch.mix(ch_meta_faa)
    combined_ch = combined_ch.mix(ch_meta_ffn)

    if (!params.skip_bakta) {
        // Run BAKTA
        bakta_out = BAKTA_BAKTA(
            ch_bakta_input,
            params.bakta_db ? ch_bakta_input.map{ _f -> file(params.bakta_db)} 
                : file("${params.beav_dir}/databases/bakta_db").exists() 
                ? ch_bakta_input.map{ _f -> file("${params.beav_dir}/databases/bakta_db")} 
                : ch_bakta_db_dl.combine(ch_bakta_input).map{ _db, _meta, _f -> file("${params.beav_dir}/databases/bakta_db")},
            [],
            [],
            [],
            []
        )
    
        // Define output channels
        ch_bakta_meta_gbff  = bakta_out.gbff
        ch_bakta_meta_gff   = bakta_out.gff
        ch_bakta_meta_fna   = bakta_out.fna
        ch_bakta_meta_faa   = bakta_out.faa
        ch_bakta_meta_ffn   = bakta_out.ffn
        ch_bakta_tiger_tuple = bakta_out.tiger
    } else {
        ch_bakta_meta_gbff  = channel.empty()
        ch_bakta_meta_gff   = channel.empty()
        ch_bakta_meta_fna   = channel.empty()
        ch_bakta_meta_faa   = channel.empty()
        ch_bakta_meta_ffn   = channel.empty()
        ch_bakta_tiger_tuple = channel.empty()  
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MAIN ELEMENTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    //
    // NF-CORE MODULE: MOBSUITE_RECON || Plasmid Characterization
    //
    if(!params.skip_mobsuite){
        mobsuiteresults = MOBSUITE_RECON(
            ch_meta_fna.mix(ch_bakta_meta_fna)
        )
        combined_ch = combined_ch.mix(mobsuiteresults.contig_report)
        combined_ch = combined_ch.mix(mobsuiteresults.biomarkers)
    }

    //
    // LOCAL MODULE: MACSYFINDER    || Secretion Systems
    //
    if(!params.skip_macsyfinder){

        macsyresults = MACSYFINDER_MACSYFINDER(
            ch_meta_faa.mix(ch_bakta_meta_faa),
            file("${params.beav_dir}/models/macsy-finder-models/TXSScan/metadata.yml").exists() 
                ? channel.fromPath("${params.beav_dir}/models/macsy-finder-models").combine(ch_filecount).map{ d, _m, _f -> d } 
                : mf_db.combine(ch_filecount).map{ db, _m, _f -> db }
        )

        combined_ch = combined_ch.mix(macsyresults.tsv)

    }

    //
    // LOCAL MODULE: DEFENSEFINDER  || Defense Systems
    //
    if(!params.skip_defensefinder){
        defenseout = DEFENSEFINDER_DEFENSEFINDER(
            ch_meta_faa.mix(ch_bakta_meta_faa),
            file("${params.beav_dir}/models/defense-finder-models/defense-finder-models/metadata.yml").exists() 
                ? channel.fromPath("${params.beav_dir}/models/defense-finder-models").combine(ch_filecount).map{ d, _m, _f -> d } 
                : df_db.combine(ch_filecount).map{ db, _m, _f -> db }
        )
        
        combined_ch = combined_ch.mix(defenseout.tsv)
        
    }

    //
    // NF-CORE MODULE: ANTISMASH    || Secondary Metabolite Cluster Prediction
    //
    if (!params.skip_antismash){

        antismash_out = ANTISMASH_ANTISMASH(
            ch_meta_gbff.mix(ch_bakta_meta_gbff),
            params.antismash_db 
                ? ch_filecount.map{ _f -> file(params.antismash_db) } 
                : file("${params.beav_dir}/databases/antismash_db").exists() ? ch_filecount.map{ _f -> file("${params.beav_dir}/databases/antismash_db") } 
                : ch_antismash_db_dl.combine(ch_filecount).map{ db, _m, _f -> db},
            []
        )

        antismash_final = ANTISMASH_GENBANKS_TO_TABLE(
            antismash_out.gbk_results
        )

        combined_ch = combined_ch.mix(antismash_final.subset)
    }

    //
    // LOCAL MODULES: OPERONMAPPER, IDENTIFY_VGRG_CLUSTERS       || Operon Prediction and vgrG Clusters
    //
    if (params.operon_email){
        // Make copy of email variable per input
        ch_email = channel.of("${params.operon_email}").combine(ch_filecount).map{ e, _m, _f -> e }

        // Join gff to fna with basename
        ch_fna_gff = channel.empty()
        ch_fna_gff = ch_fna_gff.mix(ch_meta_fna, ch_meta_gff).groupTuple()
        ch_b_fna_gff = ch_bakta_meta_gff.mix(ch_bakta_meta_fna).groupTuple()
        ch_operonmapper_input = ch_fna_gff.mix(ch_b_fna_gff)

        operonmapperout = OPERONMAPPER_SUBMITOPERONMAPPERJOB(
            ch_operonmapper_input,
            ch_email
        )

        ch_operons = operonmapperout.tsv
        ch_vgrg_input = ch_operons.mix(ch_meta_faa, ch_meta_gff).groupTuple()

        vgrgout = IDENTIFY_VGRG_CLUSTERS(
            ch_vgrg_input,
            ch_vgrg_input.map{ _f -> file("${params.beav_dir}/models/vgrgmodels/TIGR01646.1.hmm") },
            ch_vgrg_input.map{ _f -> file("${params.beav_dir}/models/vgrgmodels/TIGR03361.1.hmm") }
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
            ch_meta_faa.mix(ch_bakta_meta_faa),
            file("${params.beav_dir}/databases/PaperBLAST_db/tmp/path.aa/curated.faa").exists() 
                        ? ch_filecount.map{ _f -> file("${params.beav_dir}/databases/PaperBLAST_db") } 
                        : gap_db.combine(ch_filecount).map{ env, _m, _file -> env},
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
            ch_tDNA_leftborder  = channel.fromPath("${params.beav_dir}/models/agromodels/T-DNA_leftborder.hmm").combine(ch_filecount).map{ t, _m, _f -> t }
            ch_tDNA_rightborder = channel.fromPath("${params.beav_dir}/models/agromodels/T-DNA_rightborder.hmm").combine(ch_filecount).map{ t, _m, _f -> t }
            ch_overdrive        = channel.fromPath("${params.beav_dir}/models/agromodels/overdrive.hmm").combine(ch_filecount).map{ t, _m, _f -> t }
        }

        ch_nodbox = channel.fromPath("${params.beav_dir}/models/agromodels/nodbox.hmm").combine(ch_filecount).map{ t, _m, _f -> t }
        ch_ttsbox = channel.fromPath("${params.beav_dir}/models/agromodels/ttsbox.hmm").combine(ch_filecount).map{ t, _m, _f -> t }


        hmmerout = HMMER_NHMMER(
            ch_meta_fna.mix(ch_bakta_meta_fna),
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
        def ag = params.agro ? params.agro : params.agrobacterium
        ch_agro = channel.of(ag).combine(ch_filecount).map{ state, _m, _f -> state }

        fuzzoutm = FUZZNUC_MAIN_ELEMENTS(
            ch_meta_fna.mix(ch_bakta_meta_fna),
            ch_agro
        )
        
        combined_ch = combined_ch.mix(fuzzoutm.tsv)
    }

    //
    // NF-CORE MODULE: BLASTN           || Nucleotide Sequence Similarity Search
    //
    if (!params.skip_blastn){
        // Populate input channels
        ch_oriT = channel.fromPath("${params.beav_dir}/oriT_db.fna").combine(ch_filecount).map{ t, _m, _f -> t }
        ch_dif  = channel.fromPath("${params.beav_dir}/dif_db.fna").combine(ch_filecount).map{ t, _m, _f -> t }
        ch_blastn_false = channel.of(false).combine(ch_filecount).map{ falsity, _meta, _file -> falsity }

        blout = BLAST_BLASTN(
            ch_meta_fna.mix(ch_bakta_meta_fna),
            ch_oriT,
            ch_dif,
            [],
            [],
            ch_blastn_false
        )

        bedout = BEDTOOLS(
            blout.out
        )

        combined_ch = combined_ch.mix(bedout.oriT)
        combined_ch = combined_ch.mix(bedout.dif)
    }

    //
    // NF-CORE MODULES: SOURMASH_SKETCH, SOURMASH_COMPARE   || Genome Sketching and Comparison
    //
    if(!params.skip_sourmash){
        sourmashout = SOURMASH_SKETCH(
            ch_meta_fna.mix(ch_bakta_meta_fna).map{ meta, fasta -> tuple(meta, [fasta]) }, 
            false
        )

        signatures_ch = sourmashout.signatures

        ch_true = channel.of(true).combine(ch_filecount).map{ t, _m, _f -> t }

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
        
        ch_tiger_dir = file("${params.beav_dir}/TIGER-TIGER2.1").exists() 
                            ? channel.fromPath("${params.beav_dir}/TIGER-TIGER2.1").combine(ch_filecount).map{ db, _m, _file -> db } 
                            : t_db.combine(ch_filecount).map { _db, _m, _f -> file("${params.beav_dir}/TIGER-TIGER2.1")}
        ch_tiger_blast_db = channel.fromPath("${params.tiger_blast_db}").combine(ch_filecount).map{ db, _m, _file -> db }

        tigerout = TIGER_TIGER( 
            ch_tiger_tuple.mix(ch_bakta_tiger_tuple),
            ch_tiger_dir,
            ch_tiger_blast_db
        )

        combined_ch = combined_ch.mix(tigerout.table)
                
    } else if (!params.skip_tiger) {
        println("""
        TIGER BLAST database not specified. A database is required for TIGER ICE analysis.
        Specify with --tiger_blast_db path/to/database
        "If you do not want to run TIGER you may safely ignore this message.
        To prevent future messages from appearing, specify --skip_tiger
        """)
    }

    //
    // NF-CORE MODULE: INTEGRONFINDER   || Integron Finder
    //
    if(!params.skip_integronfinder){
        integronout = INTEGRONFINDER( ch_meta_fna.mix(ch_bakta_meta_fna) )
        combined_ch = combined_ch.mix(integronout.summary, integronout.integrons)
    }

    //
    // NF-CORE MODULE: PHISPY           || Prophage Identification
    //
    if(!params.skip_phispy){
        phispyout = PHISPY( ch_meta_gbff.mix(ch_bakta_meta_gbff) )
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
        ch_refs = channel.fromPath("${params.beav_dir}/agrobacterium_taxa/*.fna").combine(ch_filecount).map{ r, _m, _f -> tuple([id: "references", single_end: false], r) }
        fastaniout = FASTANI(
            ch_meta_fna.mix(ch_bakta_meta_fna),    // query genomes
            ch_refs,  // reference genomes
            [],
            []
        )

        combined_ch = combined_ch.mix(fastaniout.ani)

        //
        // LOCAL MODULE: IDENTIFY TI RI PLASMIDS
        //
        ch_oncgs = channel.fromPath("${params.beav_dir}/oncogenic_plasmids").combine(ch_filecount).map{ o, _m, _f -> o }
        tiriout = IDENTIFY_TI_RI_PLASMIDS(
            ch_meta_fna.mix(ch_bakta_meta_fna),
            ch_oncgs
        )

        contigs_ch  = channel.empty()
        contigs_ch  = tiriout.contiglist
        
    } else {
        contigs_ch  = channel.empty()
    }

    ch_grouped_tuple = combined_ch.groupTuple() // Group by meta id

    combined_annotations = COMBINE_ANNOTATIONS( ch_grouped_tuple )

    circos_ch = channel.empty()
    circos_ch = circos_ch.mix(combined_annotations.finalgbk, combined_annotations.finalffn, combined_annotations.finalfaa, contigs_ch).groupTuple()

    CIRCOS_CIRCOS( circos_ch )

}
