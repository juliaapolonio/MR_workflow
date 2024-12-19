process TWOSAMPLEMR {
  """
  Run TwoSampleMR on sumstats data
  """

  label 'process_medium'

  container "${ workflow.containerEngine == 'singularity' ? 'docker://juliaapolonio/julia_2smr:latest':
            'docker.io/juliaapolonio/julia_2smr:latest' }"

  input:
    tuple path(exposure), val(meta), path(outcome)
    path(reference)

  output:
    path("*pleiotropy.txt")              , emit: pleiotropy
    path("*metrics.txt")                 , emit: metrics
    path("*heterogeneity.txt")           , emit: heterogeneity
    path("*steiger.txt")                 , emit: steiger
    path("*rsquare.txt")                 , emit: rsquare
    path("*png")                         , emit: effplot
    path("*mrpresso.txt")                , emit: mrpresso
  
  when:
    task.ext.when == null || task.ext.when

  script:
  """
  run_twosamplemr.R \\
    $exposure \\
    $outcome \\
    $reference
  """
}
