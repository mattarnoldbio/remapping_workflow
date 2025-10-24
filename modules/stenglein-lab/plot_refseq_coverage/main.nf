process PLOT_REFSEQ_COVERAGE {
  label 'process_single'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.3.2"
  }     

  input:
  path (depth)
  path (R_lib_dir)

  output:
  path "*.pdf"                       , emit: pdf
  // path "collected_coverage_plot.pdf" , emit: coverage_plot

  when:
  task.ext.when == null || task.ext.when

  script:

  def args             = task.ext.args ?: ''

  """
   plot_refseq_coverage.R $depth $R_lib_dir
  """

}
