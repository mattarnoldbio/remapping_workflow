include { EXTRACT_STRAND_BIAS   } from '../../modules/stenglein-lab/extract_strand_bias'
include { BAM_TO_SAM            } from '../../modules/stenglein-lab/bam_to_sam'
include { PREPEND_TSV_WITH_ID   } from '../../modules/stenglein-lab/prepend_tsv_with_id'
                                                                                
workflow QUANTIFY_STRAND_BIAS {
                                                                                
  take:
   bam                      // [meta, bam]
   R1_antisense_orientation // boolean
 
  main:                                                                         

   BAM_TO_SAM(bam)

   // TODO: optionally make R1 orientation configurable per sample
   //       instead of all for one sample
 
   EXTRACT_STRAND_BIAS(BAM_TO_SAM.out.sam, R1_antisense_orientation)
   PREPEND_TSV_WITH_ID(EXTRACT_STRAND_BIAS.out.txt)                              

  emit:

   strand_bias = PREPEND_TSV_WITH_ID.out.prepended
 
}         
