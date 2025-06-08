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
    // Debug: check what type genetic_map is
    // println "genetic_map type: ${genetic_map.getClass()}"
    // def gmappath = file(genetic_map).toAbsolutePath().toString()
    def args                       =   task.ext.args ?: ''
    def prefix                     =   task.ext.prefix ?: "${meta.id}"
    if (!(args ==~ /.*--seed.*/)) {args += " --seed=1"}

    """
    \$(Rscript -e "d=QUILT::quilt_chunk_map('${chr}', '${genetic_map}', ${min_bp}, ${min_cm})[,3];write.table(cbind(do.call(rbind,strsplit(d, '[:|-]')), '${vcfpath}', '${indexpath}', '${gmappath}'), file='${prefix}.csv', sep=',', col.names=c('chrom', 'start', 'end', 'refpanel_vcf', 'refpanel_vcf_index', 'genetic_map'),row.names=FALSE,quote=FALSE)")
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-quilt2: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"QUILT\\")))")
    END_VERSIONS
    """
}

