#!/usr/bin/env nextflow

include { REMAPPING_WORKFLOW   } from './subworkflows/stenglein-lab/remapping_workflow'

// main named workflow
workflow MAIN_WORKFLOW {

    // main remapping workflow
    REMAPPING_WORKFLOW ()
}

// entry workflow
// https://www.nextflow.io/docs/latest/reference/syntax.html#workflow
workflow {
    MAIN_WORKFLOW ()
}

