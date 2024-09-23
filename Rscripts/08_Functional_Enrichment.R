##########################################
# Name: Priyansh Srivastava ##############
# Email: spriyansh29@gmail.com ###########
# Web: https://www.metapriyansh.com/ #####
# Script Title: Functional Enrichment ####
##########################################

## Load Custom Function
source("Rscripts/helperFunctions/go_enrichment.R")

## Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggpubr))

## Load background
count_matrix <- read.table("00_RData/tables/06_A_Count_Table.txt",
                           sep = "\t", row.names = 1, header = TRUE)
all_genes <- rownames(count_matrix) 

## Load the peak annotations
peak_metadata <- read_tsv("00_RData/tables/06_B_Peak_Metadata_Table.txt") %>% as.data.frame()

## Get all ensemble IDs
entrez_ids <- peak_metadata$ENTREZID
entrez_ids <- entrez_ids[!is.na(entrez_ids)]
background_entrez_ids <- unique(entrez_ids)

## Load Significant genes

## Significant_Data (Antibody*Treatment)
sig_df <- read.table("00_RData/tables/07_DESeq2_Tabs/antibody_sig.txt",
                     header = TRUE, sep = "\t", row.names = 1)

## UP
sig_df_up <- sig_df[sig_df$log2FoldChange >= 1.5,, drop=FALSE]
gene_list_up <- rownames(sig_df_up)
peak_metadata_up <- peak_metadata[peak_metadata$rownames %in% gene_list_up,,drop=FALSE]
entrez_ids_up <- peak_metadata_up$ENTREZID
entrez_ids_up <- entrez_ids_up[!is.na(entrez_ids_up)]
entrez_ids_up <- unique(entrez_ids_up)

## Down
sig_df_down <- sig_df[sig_df$log2FoldChange <= -1.5,, drop=FALSE]
gene_list_down <- rownames(sig_df_down)
peak_metadata_down <- peak_metadata[peak_metadata$rownames %in% gene_list_down,,drop=FALSE]
entrez_ids_down <- peak_metadata_down$ENTREZID
entrez_ids_down <- entrez_ids_down[!is.na(entrez_ids_down)]
entrez_ids_down <- unique(entrez_ids_down)

## Run GO enrichment
go_BP_up <- go_enrichment(gene_list = entrez_ids_up,
              background = background_entrez_ids,
              ont = "BP",
              nterms = 10)
go_BP_down <- go_enrichment(gene_list = entrez_ids_down,
                          background = background_entrez_ids,
                          ont = "BP",
                          nterms = 10)
go_MF_up <- go_enrichment(gene_list = entrez_ids_up,
                       background = background_entrez_ids,
                       ont = "MF",
                       nterms = 10)
go_MF_down <- go_enrichment(gene_list = entrez_ids_down,
                       background = background_entrez_ids,
                       ont = "MF",
                       nterms = 10)
go_CC_up <- go_enrichment(gene_list = entrez_ids_up,
                       background = background_entrez_ids,
                       ont = "CC",
                       nterms = 10)
go_CC_down <- go_enrichment(gene_list = entrez_ids_down,
                       background = background_entrez_ids,
                       ont = "CC",
                       nterms = 10)
antibody <- ggarrange(
  go_BP_up$dot, go_MF_up$dot, go_CC_up$dot,
  go_BP_down$dot, go_MF_down$dot, go_CC_down$dot,
  labels = c("A. Biological Process(~Antibody)(UP)", "B. Molecular Functions(~Antibody)(UP)", "C. Cellular Compartment(~Antibody)(UP)",
             "D. Biological Process(~Antibody)(DOWN)", "D. Molecular Functions(~Antibody)(DOWN)", "E. Cellular Compartment(~Antibody)(DOWN)"),
  ncol = 3, nrow = 2
)
antibody

ggsave(plot = antibody,
       filename = "00_RData/images/07_Antibody_Enrichment_Results.png",
       dpi = 600, height = 12, width = 20)


## Significant_Data (Antibody*Treatment)
sig_df <- read.table("00_RData/tables/07_DESeq2_Tabs/antibody_and_interact_treatment_sig.txt",
                     header = TRUE, sep = "\t", row.names = 1)

## UP
sig_df_up <- sig_df[sig_df$log2FoldChange >= 1.5,, drop=FALSE]
gene_list_up <- rownames(sig_df_up)
peak_metadata_up <- peak_metadata[peak_metadata$rownames %in% gene_list_up,,drop=FALSE]
entrez_ids_up <- peak_metadata_up$ENTREZID
entrez_ids_up <- entrez_ids_up[!is.na(entrez_ids_up)]
entrez_ids_up <- unique(entrez_ids_up)

## Down
sig_df_down <- sig_df[sig_df$log2FoldChange <= -1.5,, drop=FALSE]
gene_list_down <- rownames(sig_df_down)
peak_metadata_down <- peak_metadata[peak_metadata$rownames %in% gene_list_down,,drop=FALSE]
entrez_ids_down <- peak_metadata_down$ENTREZID
entrez_ids_down <- entrez_ids_down[!is.na(entrez_ids_down)]
entrez_ids_down <- unique(entrez_ids_down)

## Run GO enrichment
go_BP_up <- go_enrichment(gene_list = entrez_ids_up,
                          background = background_entrez_ids,
                          ont = "BP",
                          nterms = 10)
go_BP_down <- go_enrichment(gene_list = entrez_ids_down,
                            background = background_entrez_ids,
                            ont = "BP",
                            nterms = 10)
go_MF_up <- go_enrichment(gene_list = entrez_ids_up,
                          background = background_entrez_ids,
                          ont = "MF",
                          nterms = 10)
go_MF_down <- go_enrichment(gene_list = entrez_ids_down,
                            background = background_entrez_ids,
                            ont = "MF",
                            nterms = 10)
go_CC_up <- go_enrichment(gene_list = entrez_ids_up,
                          background = background_entrez_ids,
                          ont = "CC",
                          nterms = 10)
go_CC_down <- go_enrichment(gene_list = entrez_ids_down,
                            background = background_entrez_ids,
                            ont = "CC",
                            nterms = 10)

antibody_and_interact_treatment <- ggarrange(
  go_BP_up$dot, go_MF_up$dot, go_CC_up$dot,
  go_BP_down$dot, go_MF_down$dot, go_CC_down$dot,
  labels = c("A. Biological Process(~Antibody*Treatment)(UP)", "B. Molecular Functions(~Antibody*Treatment)(UP)", "C. Cellular Compartment(~Antibody*Treatment)(UP)",
             "D. Biological Process(~Antibody*Treatment)(DOWN)", "D. Molecular Functions(~Antibody*Treatment)(DOWN)", "E. Cellular Compartment(~Antibody*Treatment)(DOWN)"),
  ncol = 3, nrow = 2
)
antibody_and_interact_treatment

ggsave(plot = antibody_and_interact_treatment,
       filename = "00_RData/images/07_Antibody_Interact_Treatment_Enrichment_Results.png",
       dpi = 600, height = 12, width = 22)
