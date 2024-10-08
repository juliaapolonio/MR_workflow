## Introduction

**juliaapolonio/MR_workflow** is a pipeline for Mendelian Randomization and sensitivity analysis between a given GWAS sumstats and QTL data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. 
As a future improvement, when possible, the local modules will be submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Pipeline summary

```mermaid
graph TD
    A[Input] --> B[GSMR]
    B --> C[Filtering]
    C --> D[TwoSampleMR]
    C --> E[Coloc]
    E --> F[Aggregration of Results]
    D --> F
```

### Generalized Summary Mendelian Randomization (GSMR)

This is the main part of the process. It runs [GSMR](https://yanglab.westlake.edu.cn/software/gcta/#MendelianRandomisation) for all Exposures vs the Outcome and returns number of IVs, betas, SEs and p-values for each exposure.

### Significant gene calculation and filtering

With the results from GSMR, this module calculates the FDR p-value for each gene and filters by it and the number of IVs.

### Two Sample MR (2SMR)

The following two modules are run only for the FDR-significant genes of GSMR.

[Two Sample MR](https://mrcieu.github.io/TwoSampleMR/) is an R package that performs Mendelian Randomization and sensitivity analysis. The workflow is configured to run the following 2SMR tests:
- Inverse Variance Weighted regression;
- Simple Median regression;
- Simple mode regression;
- MR Egger regression;
- Heterogeneity Egger;
- Heterogeneity Inverse Variance Weighted;
- Steiger direction test;
- Pleiotropy Egger intercept.

### Coloc

[Coloc](https://chr1swallace.github.io/coloc/) is an R package to run colocalization analysis. For this workflow, the informations retrieved from Coloc are:
- H3;
- H4;
- Most probable causal variant.

### Merge all results and find candidate drug targets

An R script that collects the outputs from the previous modules and merges everything in a .csv file.

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command (WORK IN PROGRESS):

```bash
nextflow run juliaapolonio/MR_workflow -profile test,YOURPROFILE --outdir <OUTDIR>
```

Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

> - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
> - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
> - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
> - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

- Start running your own analysis!

```bash
nextflow run juliaapolonio/MR_workflow \
  --exposure samplesheet.csv \
  --outdir <OUTDIR> \
  --ref_folder 1000G_reference \
  --outcome outcome_file \
  -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
```

## Databases and references

MR_workflow needs 3 inputs to run:
- A reference folder;
- An Exposure sample sheet;
- An Outcome file.

The code to download and format reference data (in the example is European 1000Genomes) is located in `bin/setup_bfile.sh`.
Both Exposure and Outcome files should follow the [GCTA-Cojo format](https://yanglab.westlake.edu.cn/software/gcta/#COJO). The Exposure file should be separated by one gene per file.
Neither the Exposure nor the Outcome files should contain multi-allelic SNPs and the frequency is the Minor Allele Frequency (MAF).

## Outputs

If successfully run, the workflow should give two files as an output:

- `mr_merged_results.csv` should contain all analyses results for each GSMR significant gene;
- `significant_genes.txt` should give a gene list of all genes that fill the criteria.

## Credits

juliaapolonio/MR_workflow was co-authored by Julia Apolonio and João Cavalcante, under the supervision of Dr. Vasiliki Lagou.


## Citations

> **Causal associations between risk factors and common diseases inferred from GWAS summary data.**
> 
> Zhihong Zhu, Zhili Zheng, Futao Zhang, Yang Wu, Maciej Trzaskowski, Robert Maier, Matthew R. Robinson, John J. McGrath, Peter M. Visscher, Naomi R. Wray & Jian Yang
>
> _Nature Communications_ 2018 Jan 15. doi: [10.1038/s41467-017-02317-2](https://doi.org/10.1038/s41467-017-02317-2)


> **The MR-Base platform supports systematic causal inference across the human phenome.**
>
> Hemani G, Zheng J, Elsworth B, Wade KH, Baird D, Haberland V, Laurin C, Burgess S, Bowden J, Langdon R, Tan VY, Yarmolinsky J, Shihab HA, Timpson NJ, Evans DM, Relton C, Martin RM, Davey Smith G, Gaunt TR, Haycock PC, The MR-Base Collaboration.
>
> _eLife_ 2018 Jul. doi: [10.7554/eLife.34408](https://doi.org/10.7554/eLife.34408)


> **Eliciting priors and relaxing the single causal variant assumption in colocalisation analyses**
>
> Chris Wallace
>
> _PLOS Genetics_ 2020 Apr 20. doi: [10.1371/journal.pgen.1008720](https://doi.org/10.1371/journal.pgen.1008720)


> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
