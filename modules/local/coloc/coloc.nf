process COLOC {
  """
  Run Coloc on sumstats data and QTL
  """

  label 'process_medium'
  label 'ERRO'

  container "juliaapolonio/coloc:5.2.3"

  input:
    each(reads)
    path(outcome)

  output:
    path("*coloc*")        , emit: merged_coloc
    path("*png")        , emit: plots

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  colocalization.R \\
    $reads \\
    $outcome \\

  """
}
