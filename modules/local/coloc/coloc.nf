process COLOC {
  """
  Run Coloc on sumstats data and QTL
  """

  label 'process_medium'

  container "${ workflow.containerEngine == 'singularity' ? 'docker://juliaapolonio/coloc:5.2.3':
            'docker.io/juliaapolonio/coloc:5.2.3dev' }"

  input:
    tuple path(ref_bim), path(reads), val(meta), path(outcome)

  output:
    path("*_coloc.txt")      , emit: merged_coloc
    path("*_coloc.png")      , emit: plots_coloc
    path("*_regional.png")   , emit: plots_regional, optional: true

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  colocalization.R \\
    $reads \\
    $outcome \\
    $ref_bim \\

  """
}
