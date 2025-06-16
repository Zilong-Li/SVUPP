#!/usr/bin/env nextflow

include { CUTESV2 } from '../../modules/cuteSV2/force_calling/main'


workflow CUTESV2_RUN {

    take:
    samples_csv      // csv file
    ch_labels_csv    // channel: path(csv)
    
    main:

    // Parse sample sheet
    ch_samples = Channel
        .fromPath(samples_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> 
            def bam = file(row.bam, checkIfExists: true)
            def bai = file(row.bai, checkIfExists: true)
            def fasta = file(row.fasta, checkIfExists: true)
            def fai = file(row.fai, checkIfExists: true)
            tuple(row.sample, bam, bai, fasta, fai)
        }

    ch_labels = ch_labels_csv
        .splitCsv(header: true)
        .map { row -> 
            def label = file(row.label, checkIfExists: true)
            tuple(row.sample, label)
        }

    ch_input = ch_samples.join(ch_labels, by: 0)

    CUTESV2(ch_input, params.svfile)
    ch_vcf      = CUTESV2.out.vcf
    ch_tbi      = CUTESV2.out.tbi
    ch_versions = CUTESV2.out.versions

    emit:
    vcf      = ch_vcf      // channel: [ val(meta), path(vcf) ]
    tbi      = ch_tbi      // channel: [ val(meta), path(vcf) ]
    versions = ch_versions // channel: [ path(versions.yml) ]

}

