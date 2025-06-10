#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { QUILT2_PREPARE_RDATA } from './subworkflows/quilt2_prepare/main'
include { QUILT2_RUN } from './subworkflows/quilt2_run/main'
include { CUTESV2_RUN } from './subworkflows/cutesv2_run/main'


workflow {

    if (params.refpanel) {
        // Parse reference sheet
        ch_refpanel = Channel
            .fromPath(params.refpanel, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                [
                    [id: row.chrom],
                    row.chrom,
                    file(row.genetic_map, checkIfExists: true),
                    file(row.vcf, checkIfExists: true),
                    file(row.vcf_index, checkIfExists: true),
                    row.genetic_map
                ]
            }
            .unique { it[0].id }

        ch_params = Channel.of([params.minbp, params.mincm]) // params for chunking

        QUILT2_PREPARE_RDATA(ch_refpanel.combine(ch_params))
        ch_refdata = QUILT2_PREPARE_RDATA.out.csv
    } else {
        ch_refdata = Channel.fromPath(params.refdata, checkIfExists: true)
    }

    if (params.labels) {
        ch_labels = Channel.fromPath(params.labels, checkIfExists: true)
    } else {
        QUILT2_RUN(params.samples, ch_refdata)
        ch_labels =  QUILT2_RUN.out.labels
    }
    
    def samples = params.samples2 ?: params.samples
    CUTESV2_RUN(samples, ch_labels)

}

