#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { QUILT2_PREPARE_RDATA } from './subworkflows/quilt2_prepare/main'
include { QUILT2_RUN } from './subworkflows/quilt2_run/main'
include { CUTESV2_RUN } from './subworkflows/cutesv2_run/main'


workflow {

    if (params.read_labels) {
        ch_labels = Channel.fromPath(params.read_labels, checkIfExists: true)
    } else {
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
                        file(row.vcf_index, checkIfExists: true)
                    ]
                }
                .unique { it[0].id }

            ch_params = Channel.of([params.minbp, params.mincm]) // params for chunking

            QUILT2_PREPARE_RDATA(ch_refpanel.combine(ch_params))
            // save the paths in a CSV
            ch_refdata = QUILT2_PREPARE_RDATA.out.rdata
                .map { meta, rdata -> "${meta.id},${rdata}" }
                .collectFile(
                    name: 'prepared_reference_rdata.csv',
                    seed: 'chunk_id,refpanel_rdata',
                    newLine: true,
                    storeDir: params.outdir,
                    sort: false
                )
        } else {
            ch_refdata = Channel.fromPath(params.refdata, checkIfExists: true)
        }

        ch_samples = Channel.fromPath(params.samples, checkIfExists: true)
        QUILT2_RUN(ch_samples, ch_refdata)
        // save the paths in a CSV
        ch_labels = QUILT2_RUN.out.labels
            .transpose()
            .map { meta, tsv ->
                def bn = tsv.baseName.tokenize('.')
                [bn[0], tsv]
            }
            .map { id, tsv -> "${id},${tsv}" }
            .collectFile(
                name: 'samples_read_labels.csv',
                seed: 'sample,label',
                newLine: true,
                storeDir: params.outdir,
                sort: false
            )

    }
    
    def samples = params.samples2 ?: params.samples
    ch_samples = Channel.fromPath(samples, checkIfExists: true)
    CUTESV2_RUN(ch_samples, ch_labels)

}

