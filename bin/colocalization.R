#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

qtl_path <- args[1]
gwas_path <- args[2]
ref_bim <- args[3]

# Run colocalization
library(coloc)
library(dplyr)
library(ggplot2)
library(locuszoomr)
library(cowplot)

eqtl <- read.table(file=qtl_path, header=T, as.is=T); head(eqtl) 
gwas <- read.table(file=gwas_path, header=T, as.is=T); head(gwas)

eqtl$varbeta_eqtl <- eqtl$se^2
gwas$varbeta_gwas <- gwas$se^2

input <- merge(eqtl, gwas, by="SNP", all=FALSE, suffixes=c("_eqtl","gwas")); head(input)

dataset1=list(beta=input$bgwas, varbeta=input$varbeta_gwas, snp=input$SNP, type="quant", N=nrow(gwas), MAF=input$freqgwas)
dataset2=list(beta=input$b_eqtl, varbeta=input$varbeta_eqtl, snp=input$SNP, type="quant", N=nrow(eqtl), MAF=input$freq_eqtl)

result <- coloc.abf(dataset1,dataset2)


# Plot colocalization results
name <- sub("_GSMR\\.txt$", "", qtl_path)
png(filename = paste0(name,"_coloc.png"))
sensitivity(result,rule="H4 >= 0.8")
dev.off()

# Write H3 and H4 and the causal SNP in a file
causal_snp <- result[[2]][which(result[[2]][11]==max(result[[2]][11])),1]

h3 = result[[1]][5]
h4 = result[[1]][6]

line <- paste(name, h3, h4, causal_snp, sep = "\t")
write(line,file=paste0(name,"_coloc.txt"))

# Regional plot using locuszoom
if (h4 >= 0.8){

  bim <- vroom::vroom(ref_bim, col_names=c("chr", "SNP", "A", "position", "B", "C"))[,c("chr", "SNP", "position")] %>% as.data.frame()
  result_coloc <- as.data.frame(result[["results"]])
  pval <- input %>% dplyr::select("SNP", "pgwas", "p_eqtl");
  
  result_coloc <- inner_join(result_coloc, pval, by = c("snp"="SNP"));
  result_coloc <- inner_join(result_coloc, bim, by = c("snp"="SNP"));
  input_locuszoom <- subset(result_coloc, select = c("chr", "position", "snp", "pgwas", "p_eqtl"));
  causal_snp <- result[[2]][which(result[[2]][11]==max(result[[2]][11])),1]

  if (require(EnsDb.Hsapiens.v75)) {
    loc_gwas <- locus(data = input_locuszoom, gene = name, flank = 1e5,
                      ens_db = "EnsDb.Hsapiens.v75", chrom = "chr", pos = "position",
                      p = "pgwas")
    loc_gwas <- link_LD(loc_gwas, token = "5c4d1f5eeb21")
    loc_qtl <- locus(data = input_locuszoom, gene = name, flank = 1e5,
                    ens_db = "EnsDb.Hsapiens.v75", chrom = "chr", pos = "position",
                    p = "p_eqtl")
    loc_qtl <- link_LD(loc_qtl, token = "5c4d1f5eeb21")
  }

  g <- gg_genetracks(loc_gwas, highlight = name)
  pg <- gg_scatter(loc_gwas, labels = causal_snp, nudge_x = 0.2, nudge_y = 0.2) + 
    labs(title = "GWAS")
  pq <- gg_scatter(loc_qtl) +
    labs(title = "QTL")

  plot <- plot_grid(pq, pg, g, ncol = 1, rel_heights = c(2, 2, 1), align = "v")

  ggsave(paste0(name,"_regional.png"), plot, width = 3000, height = 2000, units = "px")
  
}
