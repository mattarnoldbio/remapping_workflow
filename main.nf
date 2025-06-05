#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { REMAPPING_WORKFLOW } from './subworkflows/stenglein-lab/remapping_workflow'

workflow {
  main:
    REMAPPING_WORKFLOW ()
}
