
process CUTESV2 {
    tag "$sample"
    label 'process_medium'

    conda "${moduleDir}/environment.yaml"

    input:
    tuple val(sample), path(bam), path(bai), path(fasta), path(labels)
    path(svfile)
    
    output:
    tuple val(sample), path("*.vcf.gz"),              emit: vcf
    tuple val(sample), path("*.vcf.gz.tbi"),          emit: tbi,   optional:true
    path "versions.yml",                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                        =   task.ext.args ?: ''
    def prefix                      =   task.ext.prefix ?: ''

    """
    cuteSV ${bam} ${fasta} ${sample}.vcf . \\
    --sample ${sample} \\
    --threads ${task.cpus} \\
    --min_support 1 \\
    --diff_ratio_merging_INS 0.3 \\
    --max_cluster_bias_DEL 100 \\
    --diff_ratio_merging_DEL 0.3 \\
    -read_hap1_prob ${labels} \\
    -Ivcf ${svfile} 
    $args

    bcftools sort ${sample}.vcf -o ${sample}.vcf.gz
    tabix ${sample}.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(echo \$(python --version 2>&1) | sed 's/^Python //' )
        cuteSV: \$(echo \$(cuteSV --version 2>&1) | sed 's/^cuteSV //' )
    END_VERSIONS
    """
}
