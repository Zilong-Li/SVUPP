process QUILT2_PREPARE_CHUNK {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-quilt:2.0.1--r44h503566f_1':
        'biocontainers/r-quilt:2.0.1--r44h503566f_1' }"

    input:
    tuple val(meta), val(chr), path(genetic_map), val(vcfpath), val(indexpath), val(gmappath), val(min_bp), val(min_cm)

    output:
    path("*.csv"),                         emit: csv
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                       =   task.ext.args ?: ''
    def prefix                     =   task.ext.prefix ?: "${meta.id}"
    if (!(args ==~ /.*--seed.*/)) {args += " --seed=1"}

    """

    Rscript -e "
    library(QUILT)
    
    # Generate chunk map
    chunk_map <- quilt_chunk_map('${chr}', '${genetic_map}', ${min_bp}, ${min_cm})
    
    # Extract and parse chunk information
    chunks <- chunk_map[, 3]
    chunk_parts <- do.call(rbind, strsplit(chunks, '[:|-]'))
    
    # Create output data frame
    output_df <- data.frame(
        chrom = chunk_parts[, 1],
        start = chunk_parts[, 2], 
        end = chunk_parts[, 3],
        refpanel_vcf = '${vcfpath}',
        refpanel_vcf_index = '${indexpath}',
        genetic_map = '${gmappath}',
        stringsAsFactors = FALSE
    )
    
    # Write to CSV
    write.table(
        output_df,
        file = '${prefix}.csv',
        sep = ',',
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE
    )
    "

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-quilt2: \$(Rscript -e "cat(as.character(utils::packageVersion('QUILT')))")
    END_VERSIONS
    """
}



