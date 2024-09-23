process GENE_LIST {
  """
  Move GSMR significant genes to another folder according to genelist
  """

  label 'process_medium'
  label 'ERRO'

  container "juliaapolonio/coloc:5.2.3"

  input:
    path(exposures)
    path(gene_list)

  output:
    path("filtered/*")            , emit: filtered

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  move_sign_genes.sh \\
    $gene_list
  """
}
