/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { fromSamplesheet; validateParameters } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { TWOSAMPLEMR } from "./modules/local/twosamplemr/tsmr.nf"
include { COLOC } from "./modules/local/coloc/coloc.nf"
include { GCTA_GSMR } from "./modules/local/gsmr/gsmr.nf"
include { GCTA_MERGE_ERR } from "./modules/local/gcta_merge_err/gcta_merge_err.nf"
include { RESULT } from "./modules/local/merge_results/results.nf"
include { PREPROCESS } from "./modules/local/preprocess_fstats/preprocess.nf"
include { PROCESS_REF } from "./modules/local/process_ref/process_ref.nf"
include { PARSE_2SMR } from "./modules/local/parse_2smr/parse_2smr.nf"
include { GSMR_FILTER } from "./modules/local/gsmr_filter/gsmr_filter.nf"
include { GENE_LIST } from "./modules/local/gene_list/gene_list.nf"
include { ADD_COLUMN as ADD_H_COLUMN } from "./modules/local/add_column/add_column.nf"
include { ADD_COLUMN as ADD_S_COLUMN } from "./modules/local/add_column/add_column.nf"
include { ADD_COLUMN as ADD_P_COLUMN } from "./modules/local/add_column/add_column.nf"
include { ADD_COLUMN as ADD_M_COLUMN } from "./modules/local/add_column/add_column.nf"
include { FINAL_REPORT } from "./modules/local/final_report/final_report.nf"
include { UNTAR } from "./modules/nf-core/untar/main.nf"
include { UNTAR as UNTAR_REF } from "./modules/nf-core/untar/main.nf"
include { MERGE } from "./modules/local/merge_file/merge_file.nf"
include { R_LIFT } from "./modules/local/r_lift/r_lift.nf"
include { RENDER_REPORT } from "./modules/local/render_report/main.nf"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow {

    if(params.ref){

        ref = params.ref

        Channel
            .fromPath("${ref}/*.bim") // Filter files ending with .bim
            .set { bim_files }
    }

    if(params.run_eqtlgen) {
        PREPROCESS ()
        
    }

    if(params.run_ukb) {
        reads = Channel.fromPath(params.ukb_path)
        UNTAR (
            reads
        )

        MERGE (
            UNTAR.out.untar
        )

        R_LIFT (
            MERGE.out.merged,
            sumstats
        )
    }

    if(params.create_ref) {
        PROCESS_REF (
           params.psam_url,
           params.pgen_url,
           params.pvar_url
            )
    }
    
    data = Channel.fromSamplesheet("exposures")
    outcomes = Channel.fromSamplesheet("outcomes")

    if(!params.ref && params.run_vignette){
        zenodo_ref = [[id: 'zenodo'], file(params.zenodo_link)]
        UNTAR_REF (
            zenodo_ref
        )

        ref = UNTAR_REF.out.untar.map { it[1] } + "/ref/"

        bim_files = UNTAR_REF.out.untar
            .map { meta, path -> 
                def bim = file("${path}/ref/*.bim")
                if (bim.isEmpty()) {
                    error "No .bim file found in ${path}/ref/"
                }
                [meta, bim[0]]  // We're assuming there's only one .bim file
            }
            .map { it[1] }
        
    }

     if(!params.ref && params.run_replication){
         zenodo_ref = [[id: 'zenodo'], file(params.zenodo_link)]
         UNTAR_REF (
             zenodo_ref
         )
 
         ref = UNTAR_REF.out.untar.map { it[1] } + "/ref/"
 
         bim_files = UNTAR_REF.out.untar
             .map { meta, path ->
                 def bim = file("${path}/ref/*.bim")
                 if (bim.isEmpty()) {
                     error "No .bim file found in ${path}/ref/"
                 }
                 [meta, bim[0]]  // We're assuming there's only one .bim file
             }
             .map { it[1] }
 
     }

    data.combine(outcomes).set { og_combinations }
    GCTA_GSMR (
      og_combinations,
	  ref,
    )

    GCTA_MERGE_ERR(GCTA_GSMR.out.gsmr_err.collect())

    GSMR_FILTER (
	    GCTA_GSMR.out.gsmr_res.collect()
	    )

    GENE_LIST (
            data.collect { meta, path -> path },
	    GSMR_FILTER.out.genelist
	    )

    GENE_LIST.out.filtered.flatten().combine(outcomes).set { combinations }
    TWOSAMPLEMR (
            combinations,
            ref
            )

    bim_files.combine(combinations).set { coloc_combinations }
    
    COLOC (
            coloc_combinations
	    )

    COLOC.out.merged_coloc
         .collectFile(name: 'concatenated_coloc.csv', storeDir: "${params.outdir}/collected_files/")
         .set { concatenated_coloc }

    TWOSAMPLEMR.out.mrpresso
        .collectFile(name: 'concatenated_mrpresso.txt', storeDir: "${params.outdir}/collected_files/")
        .set { concatenated_mrpresso }
    
    TWOSAMPLEMR.out.metrics
        .collectFile(name: 'concatenated_metrics.txt', storeDir: "${params.outdir}/collected_files/")
        .set { concatenated_metrics }

    TWOSAMPLEMR.out.steiger
        .collectFile(name: 'concatenated_steiger.txt', storeDir: "${params.outdir}/collected_files/")
        .set { concatenated_steiger }

    TWOSAMPLEMR.out.heterogeneity
        .collectFile(name: 'concatenated_heterogeneity.txt', storeDir: "${params.outdir}/collected_files/")
        .set { concatenated_heterogeneity }

    TWOSAMPLEMR.out.pleiotropy
        .collectFile(name: 'concatenated_pleiotropy.txt', storeDir: "${params.outdir}/collected_files/")
        .set { concatenated_pleiotropy }


    RESULT (
	    concatenated_coloc,
            GSMR_FILTER.out.filtered_genes,
            concatenated_heterogeneity,
            concatenated_steiger,
            concatenated_pleiotropy,
            concatenated_metrics,
            concatenated_mrpresso
    )

    FINAL_REPORT (
            GSMR_FILTER.out.results_gsmr,
            RESULT.out.final_results            
    )

    RENDER_REPORT (
        FINAL_REPORT.out.forest_rds,
        FINAL_REPORT.out.volcano_rds,
        RESULT.out.final_results
    )

}
