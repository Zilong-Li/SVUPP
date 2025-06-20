process QUILT2_PREPARE_REFERENCE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-quilt:2.0.2--r44h503566f_0':
        'biocontainers/r-quilt:2.0.2--r44h503566f_0' }"

    input:
    tuple val(meta), val(chr), val(regions_start), val(regions_end),  path(reference_vcf_file), path(reference_vcf_file_index), path(genetic_map_file)
    val(buffer)
    val(nGen)

    output:
    tuple val(meta), path("RData/*.RData"),        emit: rdata
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                        =   task.ext.args ?: ''
    def _prefix                     =   task.ext.prefix ?: "${meta.id}"

    """
    QUILT_prepare_reference.R \\
        --chr=$chr \\
        --regionStart=$regions_start \\
        --regionEnd=$regions_end \\
        --buffer=$buffer \\
        --nGen=$nGen \\
        --reference_vcf_file=$reference_vcf_file \\
        --genetic_map_file=${genetic_map_file} \\
        --outputdir="." \\
        --use_mspbwt=TRUE \\
        --impute_rare_common=TRUE \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-quilt2: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"QUILT\\")))")
    END_VERSIONS
    """
}

