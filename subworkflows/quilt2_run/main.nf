#!/usr/bin/env nextflow

include { QUILT2_IMPUTE } from '../../modules/quilt2/impute/main'


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
    ch_versions = QUILT2_IMPUTE.out.versions

    emit:
    vcf      = ch_vcf      // channel: [ val(meta), path(vcf) ]
    versions = ch_versions // channel: [ path(versions.yml) ]
}


