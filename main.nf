#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { QUILT2_PREPARE_RDATA } from './subworkflows/quilt2_prepare/main'




// params.input = 'samplesheet.csv'

// // Parse sample sheet
// Channel
//     .fromPath(params.input)
//     .splitCsv(header: true)
//     .map { row -> 
//         def sample = [id: row.sample, single_end: false]
//         def bams = [file(row.bam, checkIfExists: true)]
//         def bais = [file(row.bai, checkIfExists: true)]
//         tuple(sample, bams, bais)
//     }
//     .unique { it -> it[0].id }  // Remove duplicates based on meta.id
//     .set { ch_samples }



// Set default if not provided
params.outdir = params.outdir ?: 'results'
params.nGen = params.nGen ?: 100
params.buffer = params.buffer ?: 500000
params.minbp = params.minbp ?: 3000000
params.mincm = params.mincm ?: 4

workflow {

    // Parse reference sheet
    Channel
        .fromPath(params.refpanel)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.chrom]
            def gmap = file(row.genetic_map, checkIfExists: true)
            tuple(meta, row.chrom, gmap, row.vcf, row.vcf_index, row.genetic_map)
        }
        .unique { it -> it[0].id }  // Remove duplicates based on meta.id
        .set { ch_refpanel }

    Channel
        .of([params.minbp, params.mincm])
        .set { ch_params }

    QUILT2_PREPARE_RDATA(ch_refpanel.combine(ch_params))
}

