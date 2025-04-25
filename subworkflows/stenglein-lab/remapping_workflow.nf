include { PARSE_MAPPING_SAMPLESHEET   } from '../../subworkflows/stenglein-lab/parse_mapping_samplesheet'
include { MARSHAL_FASTQ               } from '../../subworkflows/stenglein-lab/marshal_fastq'
include { BOWTIE2_BUILD_ALIGN         } from '../../subworkflows/stenglein-lab/bowtie2_build_align'

workflow REMAPPING_WORKFLOW {                                                    

  // define some empty channels for keeping track of stuff
  ch_versions     = Channel.empty()

  MARSHAL_FASTQ(params.fastq_dir, params.fastq_pattern)

  PARSE_MAPPING_SAMPLESHEET(params.mapping_samplesheet)

  // pull out just sample ID (ignoring single-end vs not) 
  // so we can join just based on sample ID 
  MARSHAL_FASTQ.out.reads.map{meta, reads -> [meta.id, meta, reads]}.set{reads_ch}
  PARSE_MAPPING_SAMPLESHEET.out.sample_sheet.map{meta, index -> [meta.id, index]}.set{samplesheet_ch}

  // drop just sample ID, bring back in original meta
  reads_ch.join(samplesheet_ch).map{id, meta, reads, fasta -> [meta, reads, fasta] }.set{ch_mapping}

  // run bowtie2 and align
  def save_unaligned = false
  def sort_bam = true
  BOWTIE2_BUILD_ALIGN (ch_mapping, save_unaligned, sort_bam)

  ch_versions = ch_versions.mix ( BOWTIE2_BUILD_ALIGN.out.versions )

}

