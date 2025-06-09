
process QUILT2_HAPLOTAG {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-quilt:2.0.2--r44h503566f_0':
        'biocontainers/r-quilt:2.0.2--r44h503566f_0' }"

    input:
    tuple val(meta), path(vcfs), path(rdata)

    output:
    tuple val(meta), path("${meta.id}/*.haptag.tsv"),          emit: labels
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                        =   task.ext.args ?: ''
    def prefix                      =   task.ext.prefix ?: "${meta.id}"
    if (!(args ==~ /.*--seed.*/)) {args += " --seed=101"}

    """
    ls -1v $rdata > all_files.txt
    ls -1v $vcfs > all_vcfs.txt
    bcftools query -l ${vcfs[0]} > samples.txt 
    mkdir -p ${meta.id}

    Rscript -e "
    parse_phase <- function(rda){
    load(rda)
    lapply(final_read_labels_prob, function(l) setNames(data.frame(l), c('qname', 'prob', 'hap')))
    }
    labels <- lapply(readLines('all_files.txt'), parse_phase)
    res <- lapply(seq_along(labels[[1]]), function(i) {do.call(rbind, lapply(labels, \\`[[\\`, i))})
    names(res) <- readLines('samples.txt')
    lapply(names(res), function(n) write.table(res[[n]], file = paste0('${meta.id}/',n,'.haptag.tsv'), quote=FALSE, col.names=FALSE, row.names=F))
    "

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-quilt2: \$(Rscript -e "cat(as.character(utils::packageVersion('QUILT')))")
    END_VERSIONS
    """
}
