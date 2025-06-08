#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { QUILT2_PREPARE_CHUNK } from '../../modules/quilt2/chunk/main'
include { QUILT2_PREPARE_REFERENCE } from '../../modules/quilt2/prepare_refpanel/main'


workflow QUILT2_PREPARE_RDATA {

    take:
    ch_refpanel   // channel: [val(meta), chrom, gmap_path, vcf, vcf_index, gmap_val, minbp, mincm]

    main:
    ch_versions = Channel.empty()
    def refdata = params.refdata ?: 'prepared_reference_rdata.csv'
    def outdir = params.outdir ?: 'results'

    QUILT2_PREPARE_CHUNK(ch_refpanel)
    // ch_csv = QUILT2_PREPARE_CHUNK.out.csv
    // ch_versions = QUILT2_PREPARE_CHUNK.out.versions
    
    QUILT2_PREPARE_CHUNK.out.csv
        .splitCsv(header: true)
        .map { row ->
            def region = row.chrom + "." + row.start + "." + row.end
            def meta = [id: region]
            tuple(meta, row.chrom, row.start, row.end, row.refpanel_vcf, row.refpanel_vcf_index, row.genetic_map)
        }
        .unique { it -> it[0].id }  // Remove duplicates based on meta.id
        .set { ch_chunks }
    
    QUILT2_PREPARE_REFERENCE(ch_chunks, params.buffer, params.nGen)
    // save the paths in a CSV
    ch_csv = QUILT2_PREPARE_REFERENCE.out.rdata
        .map { meta, rdata -> "${meta.id},${rdata}" }
        .collectFile(
            name: refdata,
            seed: 'chunk_id,refpanel_rdata',
            newLine: true,
            storeDir: outdir,
            sort: false
        )
    ch_versions = QUILT2_PREPARE_REFERENCE.out.versions

    emit:
    csv      = ch_csv                                    // channel: [ path(prepared_reference_rdata.csv) ]
    versions = ch_versions                               // channel: [ path(versions.yml) ]
}

