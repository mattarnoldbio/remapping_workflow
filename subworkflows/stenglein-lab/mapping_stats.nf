// include { SAMTOOLS_STATS                         } from '../../modules/nf-core/samtools/stats/main'
// include { SAMTOOLS_COVERAGE                      } from '../../modules/nf-core/samtools/coverage/main'
// include { SAMTOOLS_DEPTH                         } from '../../modules/nf-core/samtools/depth/main'

/*
  Calculate basic mapping summary statistics from bam
 */
workflow MAPPING_STATS {

 take:
  meta_bam_fasta  // [meta, bam, fasta]
  per_base_depth  // boolean: calculate per-base coverage depth values using samtools depth?

 main:

  // define some empty channels for keeping track of stuff
  ch_versions     = Channel.empty()                                               

  // ------------------
  // samtools stats
  // ------------------

  // split up input channel for nf-core samtools modules
  // meta_bam_fasta.map{meta, bam, fasta -> [meta, bam]}.set{meta_bam}
  // meta_bam_fasta.map{meta, bam, fasta -> [meta, fasta]}.set{meta_fasta}

  SAMTOOLS_STATS(meta_bam_fasta)
  ch_versions = ch_versions.mix ( SAMTOOLS_STATS.out.versions )      

  // ------------------
  // samtools coverage
  // ------------------

  SAMTOOLS_COVERAGE(meta_bam_fasta)
  ch_versions = ch_versions.mix ( SAMTOOLS_COVERAGE.out.versions )      

  // ------------------
  // samtools depth
  // ------------------
  // this outputs per-based coverage depth values
  // only run if requested to do so

  ch_depth = Channel.empty()
  if (per_base_depth) {
    // the second null input is placeholder for a possible interval bedfile
    SAMTOOLS_DEPTH(meta_bam_fasta)
    ch_depth    = ch_depth.mix    ( SAMTOOLS_DEPTH.out.tsv )
    ch_versions = ch_versions.mix ( SAMTOOLS_DEPTH.out.versions )      
  }

 emit: 
  versions      = ch_versions
  stats         = SAMTOOLS_STATS.out.stats
  insert_sizes  = SAMTOOLS_STATS.out.insert_sizes
  coverage      = SAMTOOLS_COVERAGE.out.coverage
  depth         = ch_depth

}


// from: https://github.com/nf-core/modules/blob/master/modules/nf-core/samtools/stats/main.nf
process SAMTOOLS_STATS {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.18--h50ea8bc_1' :
        'biocontainers/samtools:1.18--h50ea8bc_1' }"

    input:
    tuple val(meta), path(input), path(fasta)

    output:
    tuple val(meta), path("*.stats"), emit: stats
    tuple val(meta), path("*.insert_sizes.tsv"), emit: insert_sizes
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""

    """
    samtools \\
        stats \\
        --threads !{task.cpus} \\
        ${reference} \\
        ${input} \\
        > ${prefix}.stats

    # pull out insert sizes from mapping stats
    grep ^IS ${prefix}.stats | cut -f 2- | awk '{print "${meta.id}" "\t" \$0}' > ${prefix}.insert_sizes.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process SAMTOOLS_COVERAGE {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.18--h50ea8bc_1' :
        'biocontainers/samtools:1.18--h50ea8bc_1' }"

    input:
    tuple val(meta), path(input), path(fasta)

    output:
    tuple val(meta), path("*.txt"), emit: coverage
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        coverage \\
        $args \\
        -o ${prefix}.coverage.txt \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """
}

process SAMTOOLS_DEPTH {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.18--h50ea8bc_1' :
        'biocontainers/samtools:1.18--h50ea8bc_1' }"

    input:
    tuple val(meta), path(bam), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        depth \\
        --threads ${task.cpus-1} \\
        $args \\
        -o ${prefix}.depth.tsv \\
        $bam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}


