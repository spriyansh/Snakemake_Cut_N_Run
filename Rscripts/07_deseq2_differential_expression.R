#######################################
# Name: Priyansh Srivastava ##############
# Email: spriyansh29@gmail.com ###########
# Web: https://www.metapriyansh.com/ #####
# Script Title: Differential Expression ##
##########################################

## Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(DESeq2))

## Load Custom FUnction
source("Rscripts/helperFunctions/runDESeq2.R")

## Set the I/O
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
outTab_dir <- args[2]
outImg_dir <- args[3]

## Load the counts
count_table_raw <- as.data.frame(read.table(
  file = paste(input_dir,
    "06_A_Count_Table.txt",
    sep = "/"
  ),
  header = TRUE, sep = "\t", row.names = 1
))

## Load Peak Metadata
peak_metadata_raw <- read_tsv(
  file = paste(input_dir,
    "06_B_Peak_Metadata_Table.txt",
    sep = "/"
  ),
  col_names = TRUE
) %>% as.data.frame()

# Subset
select <- which(rowSums(count_table_raw) > 100) ## remove low count genes
count_table <- count_table_raw[select, , drop = FALSE]
peak_metadata <- peak_metadata_raw[peak_metadata_raw$rownames %in% rownames(count_table), , drop = FALSE]

# Load/Create Metadata
sample_metadata <- data.frame(
  sample_id = colnames(count_table),
  row.names = colnames(count_table)
)

# Add Antibody
sample_metadata$antibody <- unlist(lapply(str_split(sample_metadata$sample_id, pattern = "_"), function(x) {
  x[[2]]
}))
sample_metadata$time <- str_remove(unlist(lapply(str_split(sample_metadata$sample_id, pattern = "_"), function(x) {
  x[[1]]
})), "S")
sample_metadata$treatment <- str_remove(unlist(lapply(str_split(sample_metadata$sample_id, pattern = "_"), function(x) {
  x[[1]]
})), "S")

# Update
sample_metadata$time <- ifelse(sample_metadata$time %in% c("2D", "4D"), str_remove(sample_metadata$time, "D"), 0)
sample_metadata$treatment <- ifelse(sample_metadata$treatment %in% c("2D", "4D"), "treated", "untreated")

# Convert all string as factors
sample_metadata$antibody <- as.factor(sample_metadata$antibody)
sample_metadata$time <- as.factor(sample_metadata$time)
sample_metadata$treatment <- as.factor(sample_metadata$treatment)

# Simple Antibody
runDESeq2(
  count_table = count_table,
  metadata = sample_metadata,
  formula = "~antibody",
  img_out_path = outImg_dir,
  tab_out_path = outTab_dir,
  return_table = F
)

# Simple Time
runDESeq2(
  count_table = count_table,
  metadata = sample_metadata,
  formula = "~time",
  img_out_path = outImg_dir,
  tab_out_path = outTab_dir,
  return_table = F
)


# Simple Tretment
runDESeq2(
  count_table = count_table,
  metadata = sample_metadata,
  formula = "~treatment",
  img_out_path = outImg_dir,
  tab_out_path = outTab_dir,
  return_table = F
)

# Simple Antiboy+Time
runDESeq2(
  count_table = count_table,
  metadata = sample_metadata,
  formula = "~antibody+treatment",
  img_out_path = outImg_dir,
  tab_out_path = outTab_dir,
  return_table = F
)

# Simple Antiboy+Time
runDESeq2(
  count_table = count_table,
  metadata = sample_metadata,
  formula = "~antibody+time",
  img_out_path = outImg_dir,
  tab_out_path = outTab_dir,
  return_table = F
)

# Simple Antiboy+Time
runDESeq2(
  count_table = count_table,
  metadata = sample_metadata,
  formula = "~antibody*treatment",
  img_out_path = outImg_dir,
  tab_out_path = outTab_dir,
  return_table = F
)

# Write Sample Metadata
write.table(sample_metadata,
  file = paste(outTab_dir, "07_Sample_Metadata.txt", sep = "/"),
  quote = FALSE, sep = "\t", row.names = FALSE
)
