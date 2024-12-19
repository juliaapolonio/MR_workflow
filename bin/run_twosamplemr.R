#!/usr/bin/env Rscript
library(TwoSampleMR)
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

exposure_path <- args[1]
prefix_exp <- sub(".*/|_.*", "", exposure_path)
outcome_path <- args[2]
prefix_outcome <- sub(".*/|_.*", "", outcome_path)
ref <- args[3]

exp <-
  read_exposure_data(
    exposure_path,
    sep = "\t",
    snp_col = "SNP",
    beta_col = "b",
    se_col = "se",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "p",
    eaf_col = "freq",
    samplesize_col = "N"
  )

exp[, "exposure"] <- prefix_exp

exp_filtered <- exp[which(exp$pval.exposure < 0.00001), ]

mic_exp <- ieugwasr::ld_clump(
  exp_filtered |> dplyr::select(
    rsid = SNP,
    pval = pval.exposure,
    id = id.exposure
  ),
  clump_kb = 1000,
  clump_p = 5e-8,
  clump_r2 = 0.05,
  plink_bin = "/usr/local/bin/plink",
  bfile = paste0(ref, "/1KG_phase3_EUR")
) |>
  dplyr::select(-c(pval, id)) |>
  dplyr::left_join(
    exp,
    by = c("rsid" = "SNP")
  ) |>
  dplyr::rename(SNP = "rsid")

outcome <-
  read_outcome_data(
    outcome_path,
    sep = "\t",
    snp_col = "SNP",
    beta_col = "b",
    se_col = "se",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "p",
    eaf_col = "freq",
    samplesize_col = "N"
  )

outcome[, "outcome"] <- prefix_outcome


# Create harmonized dataset to run analysis
dat <- harmonise_data(exposure_dat = mic_exp, outcome_dat = outcome)


# Calculate rsquare for each association
df <- add_rsq(dat)

write.table(
  df,
  file = paste0(prefix_exp, "_rsquare.txt"),
  sep = "\t",
  row.names = F,
  quote = F
)

# Run MR regressions and remove spaces from method name
mr <- mr(dat)
mr$method <- gsub("\\s+", "_", mr$method)
write.table(
  mr[,c("exposure", "method", "nsnp", "b", "se", "pval")],
  file = paste0(prefix_exp, "_metrics.txt"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)

# Run heterogeneity test and remove spaces from method name
hetero <- mr_heterogeneity(dat)
hetero$method <- gsub("\\s+", "_", hetero$method)
write.table(
  hetero[,c("exposure", "method", "Q", "Q_df", "Q_pval")],
  file = paste0(prefix_exp, "_heterogeneity.txt"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)

# Run Steiger test
steiger <- directionality_test(dat)
write.table(
  steiger[,c("exposure", "snp_r2.exposure", "snp_r2.outcome", "correct_causal_direction", "steiger_pval")],
  file = paste0(prefix_exp, "_steiger.txt"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)

# Run pleiotropy test
pleiotropy <- mr_pleiotropy_test(dat)
write.table(
  pleiotropy[,c("exposure", "egger_intercept", "se", "pval")],
  file = paste0(prefix_exp, "_pleiotropy.txt"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)

# Run MR-PRESSO and get global p-value
mrpresso <- run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)
mrpresso_pval <- mrpresso[[1]][["MR-PRESSO results"]]["Global Test"]$`Global Test`$Pvalue
write(c(prefix_exp, mrpresso_pval), ncolumns = 2, file = paste0(prefix_exp, "_mrpresso.txt"))

# Effect plot
p1 <- mr_scatter_plot(mr, dat)

png(paste0(prefix_exp, "_effect.png"), width = 800, height = 400)
p1[[1]] + theme_classic(base_size = 12) + theme(legend.text = element_text(size = 14))
dev.off()
