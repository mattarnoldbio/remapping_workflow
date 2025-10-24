/*
 * This process splits a bam file created by mapping to multiple reference sequences 
 * into one bam file per refseq
 */

workflow SPLIT_BAM_BY_REFSEQ {

 take:
  bam_fasta       // [meta, bam, fasta]

 main:

  // define some empty channels for keeping track of stuff
  ch_versions     = Channel.empty()

  // split up input channel into separate bam and fasta channels


  // split_fasta_ch = bam_fasta.map{meta, bam, fasta -> fasta}
  // split_fasta_ch = bam_fasta

  split_fasta_ch = bam_fasta.map{meta, bam, fasta -> [meta, fasta]}
    .splitFasta( record: [id: true] )

  split_bam_fasta_ch = bam_fasta.map{meta, bam, fasta -> [meta, bam]}
    .combine(split_fasta_ch, by: 0)

  SPLIT_BAM_BY_ONE_REFSEQ(split_bam_fasta_ch)

 emit:

  per_refseq_bam = SPLIT_BAM_BY_ONE_REFSEQ.out.per_refseq_bam

}

process SPLIT_BAM_BY_ONE_REFSEQ {
   tag "$meta.id"
   label 'process_low'

   conda "bioconda::samtools=1.16.1"
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
       'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
       'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

   input:
   tuple val(meta), path(bam), val(refseq_id)

   output:
   tuple val(meta), path("*.bam", includeInputs: false), val(refseq_id),  emit: per_refseq_bam, optional: true
   path  "versions.yml",            emit: versions

   when:
   task.ext.when == null || task.ext.when

   script:
   def new_bam_name = bam.name.replaceAll(/.bam$/, ".${refseq_id.id}.bam")
   """
   # first have to sort bam
   samtools index $bam

   # pull out refseq of interest (a region in samtools parlance)
   samtools \\
       view \\
       --threads ${task.cpus-1} \\
       -o ${new_bam_name} \\
       $bam \\
       ${refseq_id.id}

   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
       samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
   END_VERSIONS
   """
}


