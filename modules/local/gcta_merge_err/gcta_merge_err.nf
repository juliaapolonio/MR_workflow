process GCTA_MERGE_ERR {
    publishDir "${params.outdir}/GCTA_GSMR/", mode: 'copy'

    input:
    path err_files

    output:
    path "gcta_error_genes.txt"

    script:
    """
    cat ${err_files} > gcta_error_genes.txt
    """
}