#!/usr/bin/env nextflow

include { GENERATE_GENBANK_FILES } from '../../../modules/local/fileparser/generategenbankfiles/main.nf'

workflow fileParser{

    take:
    ch_input

    main:
    // Check if ch_input contains a directory of files rather than one file, or if it contains a comma-separated list
    if (params.input.contains(',')) {
        input_paths = params.input.split(',').collect{f -> f.trim()}
    } else {
        input_paths = [ params.input ]
    }

    ch_input = channel.from(input_paths).map { f -> file(f) }
    
    ch_files = ch_input.flatMap { input ->
        def input_file = input.toFile()
        if (input_file.isDirectory()) {
            return input_file.listFiles().findAll { f -> f.isFile() }.collect { f -> file(f) }
        } else {
            return [ input ]
        }
    }

    // Check file types and split them into channels
    ch_filetype = ch_files.map { f ->
        def filename = f.name
        def type = 'unknown'

        if (filename.endsWith('.fasta') || filename.endsWith('.fa') || filename.endsWith('.fna')) {
            type = 'fasta'

        } else if (filename.endsWith('.gbff')) {
            type = 'gbff'

        } else if (filename.endsWith('.gbk')) {
            type = 'gbk'

        } else if (filename.endsWith('.gff3') || filename.endsWith('.gff')) {
            type = 'gff'

        } else if (filename.endsWith('.faa')) {
            type = 'faa'

        } else if (filename.endsWith('.ffn')) {
            type = 'ffn'

        } else if (filename.endsWith('.tsv')) {
            type = 'csv'

        } else {
            type = 'unknown'

        }

        tuple(f, type)
    }

    // Split channel by filetype and find files from samplesheet if present
    ch_csv     = ch_filetype
                    .filter { _f, type -> type == "csv" }
                    .map { f, __ -> f }

    ch_samplesheet_items = ch_csv.splitCsv().map{ row ->
        def meta = [id: row.sample, single_end: row.fastq_2 ? false : true]
        def fq1 = file(row.fastq_1, checkIfExists: true)
        def fq2 = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : null
        def fastqs = meta.single_end ? [fq1] : [fq1, fq2]

        return [ meta, fastqs ]
    }

    ch_fasta   = ch_filetype
                    .filter { _f, type -> type == "fasta" }
                    .map { f, __ -> f }
    ch_gbk     = ch_filetype
                    .filter { _f, type -> type == "gbk" }
                    .map { f, __ -> f }
    ch_gbff    = ch_filetype
                    .filter { _f, type -> type == "gbff" }
                    .map { f, __ -> f }
    ch_gff     = ch_filetype
                    .filter { _f, type -> type == "gff" }
                    .map { f, __ -> f }
    ch_ffn     = ch_filetype
                    .filter { _f, type -> type == "ffn" }
                    .map { f, __ -> f }
    

    // Remove gff/gbff files if they are missing a sequence
    // updated 10-5-26, now does so without a process using a buffered reader to prevent excessive memory use
    ch_complete_gffs = ch_gbff.mix(ch_gff).map{ f ->
        def has_seq = f.readLines().any { line ->
            line.contains('ORIGIN') || line.contains('##FASTA')
        }

        if (has_seq) {
            f
        } else {
            null
        }

    }.filter { it -> it != null}

    // Check if BAKTA has already run on certain files
    ch_all_inputs = channel.empty()
    ch_all_inputs = ch_all_inputs.mix(ch_fasta, ch_gbk, ch_complete_gffs, ch_ffn).map{ f -> tuple( [id: f.baseName, single_end: false], f ) }
    ch_all_inputs.mix(ch_samplesheet_items)
    ch_to_generate = ch_all_inputs.map{ f ->
        tuple([id: f.baseName, single_end: false], f)
    }

    // updated 10-5-26, removed process, use read buffer
    checked_ch = ch_all_inputs.map{ meta, f ->
        def bakta_annotated = f.readLines().any { line ->
            line.contains('Bakta')
        }

        if (bakta_annotated) {
            tuple(meta, f, true)
        } else {
            tuple(meta, f, false)
        }
    }

    ch_to_generate = checked_ch 
                        .map{ meta, f, state -> 
                            if (state == false) {
                                log.warn("File potentially already annotated with Bakta: ${f}. Skipping Bakta annotation for this file. If this is incorrect, please rerun the pipeline with the parameter --force_bakta_annotation")
                            }
                            tuple(meta, f, state)
                        }
    
    generatedfiles      = GENERATE_GENBANK_FILES(ch_to_generate)

    if (!params.force_bakta_annotation){
        ch_meta_needsbakta  = generatedfiles.fnaB
        ch_meta_genbank = generatedfiles.gbff.mix(generatedfiles.gbk)
        ch_meta_fasta  = generatedfiles.fna
        ch_meta_faa    = generatedfiles.faa
        ch_meta_ffn    = generatedfiles.ffn
        ch_meta_gff    = generatedfiles.gff
    } else {
        ch_meta_needsbakta = ch_all_inputs
    }

    ch_filecount = ch_meta_needsbakta.mix(ch_meta_genbank)
    
    emit:
    ch_meta_needsbakta  // fna, fa, or fastas that need bakta annotation
    ch_meta_fasta       // generated files that don't need bakta annotation
    ch_meta_faa         
    ch_meta_ffn
    ch_meta_gff
    ch_meta_genbank     // gbk and gbff files
    ch_filecount        // helper channel to help create instance of static file per unique input file

}