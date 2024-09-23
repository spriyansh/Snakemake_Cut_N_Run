###############################################
# Name: Priyansh Srivastava ###################
# Email: spriyansh29@gmail.com ################
# Web: https://www.metapriyansh.com/ ##########
# Script Title: Replicate Reproducibility #####
###############################################

## Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggcorrplot))

## Set the I/O
# inPath_hg38 <- "/home/priyansh/12_bin500_beds/"
# outPathImg <- "/home/priyansh/00_RData/images/"
# outPathTab <- "/home/priyansh/00_RData/tables/"

## Set the I/O
args <- commandArgs(trailingOnly = TRUE)
inPath_hg38 <- args[1]
outPathImg <- args[2]
outPathTab <- args[3]

## load the log files of samtools
frag.list <- lapply(list.files(path = inPath_hg38, pattern = "*.bed"),
  FUN = function(file) {
    return(file)
  }
)
names(frag.list) <- str_remove(frag.list, pattern = ".bed")

# Invoke null
reprod <- c()
fragCount <- NULL


for (hist in names(frag.list)) {
  if (is.null(fragCount)) {
    fragCount <- read.table(paste(inPath_hg38, frag.list[[hist]], sep = "/"), header = FALSE)
    colnames(fragCount) <- c("chrom", "bin", hist)
  } else {
    fragCountTmp <- read.table(paste(inPath_hg38, frag.list[[hist]], sep = "/"), header = FALSE)
    colnames(fragCountTmp) <- c("chrom", "bin", hist)
    fragCount <- full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
  }
}

## Compute Correlation
M <- cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs")

## Plot corrleation
corplot <- ggcorrplot(M, method = "circle", lab = TRUE)

# Save
ggsave(
  plot = corplot,
  filename = "03_A_Corrplot.png",
  path = outPathImg,
  dpi = 800, width = 12, height = 12
)

# Write tables
write.table(M,
  file = paste0(outPathTab, "03_Correlation_Matrix.txt"),
  quote = F, sep = "\t", row.names = T, col.names = T
)
