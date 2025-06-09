#!/usr/bin/env nextflow

include { QUILT2_IMPUTE } from '../../modules/quilt2/impute/main'
include { QUILT2_HAPLOTAG } from '../../modules/quilt2/haplotag/main'


workflow QUILT2_RUN {

    take:
    csv_samples   // channel: path(csv)
    csv_chunks    // channel: path(csv)
    
    main:

    // Parse the chunks
    ch_chunks = Channel
        .fromPath(csv_chunks)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.chunk_id, row.refpanel_rdata)
        }
        .unique { it[0] } 

    
    // Parse sample sheet
    ch_samples = Channel
        .fromPath(csv_samples)
        .splitCsv(header: true)
        .map { row -> 
            def batch = [id: row.batch]
            def bam = file(row.bam, checkIfExists: true)
            def bai = file(row.bai, checkIfExists: true)
            def fasta = row.fasta ? file(row.fasta, checkIfExists: true) : ""
            tuple(batch, bam, bai, fasta)
        }
        .groupTuple(by: 0)

    ch_input = ch_samples.combine(ch_chunks)

    QUILT2_IMPUTE(ch_input, params.nGen, params.buffer)
    ch_vcf      = QUILT2_IMPUTE.out.vcf
    ch_rda      = QUILT2_IMPUTE.out.rdata
    ch_versions = QUILT2_IMPUTE.out.versions

    ch_impute = ch_vcf
        .groupTuple(by: 0)
        .join(ch_rda.groupTuple(by: 0), by: 0)  // Join by sample ID (more explicit than combine)
        .filter { meta, vcf_list, rda_list -> 
            vcf_list.size() > 0 && rda_list.size() > 0 
        }
    
    QUILT2_HAPLOTAG(ch_impute)
    ch_labels    = QUILT2_HAPLOTAG.out.labels

    emit:
    vcf      = ch_vcf      // channel: [ val(meta), path(vcf) ]
    labels   = ch_labels   // channel: [ val(meta), path(ch_labels) ]
    versions = ch_versions // channel: [ path(versions.yml) ]

}

