#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(ggplot2)
library(ggrepel)
library(readr)
library(dplyr)
library(vroom)

######################## VOLCANO PLOT ###############################

gsmr_result_path <- args[1]
merged_result_path <- args[2]

#gsmr_result_path <- "../MR_workflow/results/results_gsmr.txt"
#merged_result_path <- "../MR_workflow/results/final_results.csv"

# Candidates file
is_candidate <- vroom(merged_result_path) # Output from merge module

is_candidate <- is_candidate[is_candidate$`is_candidate <- ifelse(gene %in% is_candidate$gene, TRUE, FALSE)` == TRUE,]

# Full file
mr_result <- vroom(gsmr_result_path, col_names = c("Exposure", "Outcome", "bxy", "se", "p", "nsnp", "heidi")) # Output from GSMR module before filter

mr_result <- mr_result %>%
  mutate(Exposure=gsub(x=Exposure,pattern="_GSMR", replacement="")) %>%
  mutate(is_candidate = ifelse(Exposure %in% is_candidate$gene, T, F))

mr_result <- mr_result %>% mutate("p"=as.numeric(p)) %>% filter(!is.na(p)) %>% filter(nsnp>3) %>% mutate(Q=p.adjust(p,method="BH"))

# Make plot
png("volcano.png", width = 750, height = 500)

ggplot(mr_result, mapping=aes(x=bxy, y=-log10(Q), label = Exposure))+
  geom_point(aes(size=nsnp), colour="gray")+
  geom_point(data=mr_result[mr_result$is_candidate == T,], aes(x=bxy, y=-log10(Q), size=nsnp), colour="black")+
  #geom_hline(yintercept=-log10(0.05/16000), linetype="dashed", color="red")+
  geom_text_repel(data=mr_result[mr_result$is_candidate == T,], aes(label=Exposure, size=20))+
  xlim(-0.05, 0.05)+
  theme_classic()

dev.off()

######################## FOREST PLOT ###############################

# Load and filter by number of SNPs
eqtlgen_all=vroom(gsmr_result_path, col_names = c("Exposure", "Outcome", "bxy", "se", "p", "nsnp", "heidi"))
eqtlgen_result=eqtlgen_all %>% mutate("p"=as.numeric(p)) %>% filter(!is.na(p))
eqtlgen_result=eqtlgen_result %>% filter(nsnp>3)

#FDR adjusted p values for remaining probes
eqtlgen_result=eqtlgen_result %>% mutate(Q=p.adjust(p,method="BH"))
eqtlgen_result=eqtlgen_result %>% filter(Q<0.05)

colnames(eqtlgen_result) <- c("p1", "p2", "estimate", "stderr", "p", "NSNPs")

eqtlgen_result$p1 <- stringr::str_remove(eqtlgen_result$p1, "\\_GSMR")

eqtlgen_result$upper <- eqtlgen_result$estimate+eqtlgen_result$stderr
eqtlgen_result$lower <- eqtlgen_result$estimate-eqtlgen_result$stderr

ci.lb = eqtlgen_result$lower
ci.ub = eqtlgen_result$upper

# taken from forest.default source code
level = 95
alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
vi <- ((ci.ub - ci.lb)/(2 * qnorm(alpha/2, lower.tail = FALSE)))^2
wi <- 1/sqrt(vi)
psize <- wi/sum(wi, na.rm = TRUE)
psize <- (psize - min(psize, na.rm = TRUE))/(max(psize, 
                                                 na.rm = TRUE) - min(psize, na.rm = TRUE))
eqtlgen_result$psize <- (psize * 1) + 0.5


# p-value filter to display results
sign = eqtlgen_result[eqtlgen_result$p<0.0005,]

ggplot(data=sign, aes(y=p1, x=estimate, xmin=estimate-stderr, xmax=estimate+stderr)) +
  geom_point(aes(size=psize)) +
  scale_size_continuous(range = c(0,2)) +
  geom_errorbarh(height=0.5) +
  labs(title='GSMR Analysis', x='Effect Size', y = 'Gene') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  guides(size=FALSE) +
  theme_classic()

ggsave("forest.png",units="in",width=5, height=15, dpi=300)
