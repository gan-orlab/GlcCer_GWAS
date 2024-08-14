# LocusZoom colocalization qq plot
## Set up environment
library(locuscomparer)
library(data.table)
library(cowplot)
library(ggplot2)

setwd("Desktop/glccer_temp/coloc")

## Load in summary stats
gwas_fn = fread("Nalls_coloc_input.txt")
colnames(gwas_fn) <- c("rsid", "pval")
gwas_fn$pval <- as.numeric(gwas_fn$pval)

glccer_fn = fread("GlcCer_Total_coloc_input.txt")
colnames(glccer_fn) <- c("rsid", "pval")
glccer_fn$pval <- as.numeric(glccer_fn$pval)


## Locuscompare
### typically this is done with GWAS stats and eQTL stats
### here we're doing it with GWAS stats and GlcCer stats
### note: lead snp from GlcCer not present in Nalls
tiff("Total_GlcCer_coloc.tiff", width = 2500, height = 1700, res = 200)
locuscompare(in_fn1 = gwas_fn, in_fn2 = glccer_fn, title = 'PD GWAS', title2 = 'Total GlcCer GWAS', snp = "rs11729411")
dev.off()

### Prep data
C16_fn = fread("GlcCer_C16_coloc_input.txt")
colnames(C16_fn) <- c("rsid", "pval")
C16_fn$pval <- as.numeric(C16_fn$pval)

C18_fn = fread("GlcCer_C18_coloc_input.txt")
colnames(C18_fn) <- c("rsid", "pval")
C18_fn$pval <- as.numeric(C18_fn$pval)

C20_fn = fread("GlcCer_C20_coloc_input.txt")
colnames(C20_fn) <- c("rsid", "pval")
C20_fn$pval <- as.numeric(C20_fn$pval)

C22_fn = fread("GlcCer_C22_coloc_input.txt")
colnames(C22_fn) <- c("rsid", "pval")
C22_fn$pval <- as.numeric(C22_fn$pval)

C23_fn = fread("GlcCer_C23_coloc_input.txt")
colnames(C23_fn) <- c("rsid", "pval")
C23_fn$pval <- as.numeric(C23_fn$pval)

C24_1_fn = fread("GlcCer_C24-1_coloc_input.txt")
colnames(C24_1_fn) <- c("rsid", "pval")
C24_1_fn$pval <- as.numeric(C24_1_fn$pval)

C24_fn = fread("GlcCer_C24_coloc_input.txt")
colnames(C24_fn) <- c("rsid", "pval")
C24_fn$pval <- as.numeric(C24_fn$pval)

### Save plots
C16_plot <- locuscompare(in_fn1 = gwas_fn, in_fn2 = C16_fn, title = 'PD GWAS', title2 = 'C16 GlcCer GWAS', snp = "rs11729411")
C18_plot <- locuscompare(in_fn1 = gwas_fn, in_fn2 = C18_fn, title = 'PD GWAS', title2 = 'C18 GlcCer GWAS', snp = "rs11729411")
C20_plot <- locuscompare(in_fn1 = gwas_fn, in_fn2 = C20_fn, title = 'PD GWAS', title2 = 'C20 GlcCer GWAS', snp = "rs11729411")
C22_plot <- locuscompare(in_fn1 = gwas_fn, in_fn2 = C22_fn, title = 'PD GWAS', title2 = 'C22 GlcCer GWAS', snp = "rs11729411")
C23_plot <- locuscompare(in_fn1 = gwas_fn, in_fn2 = C23_fn, title = 'PD GWAS', title2 = 'C23 GlcCer GWAS', snp = "rs11729411")
C24_1_plot <- locuscompare(in_fn1 = gwas_fn, in_fn2 = C24_1_fn, title = 'PD GWAS', title2 = 'C24:1 GlcCer GWAS', snp = "rs11729411")
C24_plot <- locuscompare(in_fn1 = gwas_fn, in_fn2 = C24_fn, title = 'PD GWAS', title2 = 'C24 GlcCer GWAS', snp = "rs11729411")

### Combine into cowplot
tiff(file="GlcCer_all_coloc_cowplot.tiff", height = 3600, width = 3200, res = 200)
cowplot::plot_grid(C16_plot + theme(axis.title.x = element_blank() ), 
                   C18_plot + theme(axis.title.x = element_blank() ), 
                   C20_plot + theme(axis.title.x = element_blank() ),
                   C22_plot + theme(axis.title.x = element_blank() ),
                   C23_plot + theme(axis.title.x = element_blank() ),
                   C24_1_plot + theme(axis.title.x = element_blank() ),
                   C24_plot,
                   nrow = 4,
                   labels = "auto")
dev.off()

## Coloc statistical analysis
library(coloc)
library(dplyr)

gc <- fread("GlcCer_Total_coloc_input2.txt")
colnames(gc) <- c("rsid","beta","se","maf","p","N")
pd <- fread("Nalls_coloc_input2.txt")
colnames(pd) <- c("rsid","beta","se","maf","p","N")

### Make varBetas
gc$varbeta = gc$se*gc$se
pd$varbeta = pd$se*pd$se

### Merge datasets
merged <- merge(gc, pd, by="rsid")
colnames(merged) <- c("rsid","beta_gc","se_gc","maf_gc","p_gc","N_gc","varbeta_gc","beta_pd","se_pd","maf_pd","p_pd","N_pd","varbeta_pd")

### Run coloc 
results <- coloc.abf(dataset1=list(beta=merged$beta_gc, varbeta=merged$varbeta_gc, p=merged$p_gc, snp=merged$rsid, MAF=merged$maf_gc, N=merged$N_gc, type="quant"), 
                     dataset2=list(beta=merged$beta_pd, varbeta=merged$varbeta_pd, p=merged$p_pd, maf=merged$maf_pd, snp=merged$rsid, type="cc"))
