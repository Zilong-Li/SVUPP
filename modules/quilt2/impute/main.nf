
process QUILT2_IMPUTE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-quilt:2.0.2--r44h503566f_0':
        'biocontainers/r-quilt:2.0.2--r44h503566f_0' }"

    input:
    tuple val(meta), path(bams), path(bais), val(fastas), val(chunkid), path(refdata)
    val(nGen)
    val(buffer)

    output:
    tuple val(meta), path("${meta.id}/*.vcf.gz"),              emit: vcf
    tuple val(meta), path("${meta.id}/*.vcf.gz.tbi"),          emit: tbi,   optional:true
    tuple val(meta), path("${meta.id}/RData/*.RData"),         emit: rdata
    tuple val(meta), path("${meta.id}/plots/*"),               emit: plots, optional:true
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                        =   task.ext.args ?: ''
    def prefix                      =   task.ext.prefix ?: "${meta.id}"
    def (chr, start, end)           =   chunkid.split('\\.')

    def fasta                       =   fastas.flatten().unique()
    def extensions                  =   bams.collect { it.extension }
    def extension                   =   extensions.flatten().unique()
    def list_command                =   extension == ["bam"]  ? "--bamlist="                       :
                                        extension == ["cram"] ? "--reference=${fasta[0]} --cramlist=" : ""
    if (!(args ==~ /.*--seed.*/)) {args += " --seed=101"}

    """
    printf "%s\\n" $bams | tr -d '[],' > all_files.txt
    mkdir -p ${meta.id}/RData
    mkdir -p ${meta.id}/plots
    QUILT2.R \\
        ${list_command}all_files.txt \\
        --chr=$chr \\
        --regionStart=$start \\
        --regionEnd=$end \\
        --nGen=$nGen \\
        --buffer=$buffer \\
        --nCores=$task.cpus \\
        --outputdir=${meta.id} \\
        --prepared_reference_filename=$refdata \\
        --output_read_label_prob=TRUE \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-quilt2: \$(Rscript -e "cat(as.character(utils::packageVersion('QUILT')))")
    END_VERSIONS
    """
}
