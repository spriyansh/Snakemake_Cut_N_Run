##################################
## Author: Priyansh Srivastava ###
## Email: spriyansh29@gmail.com ##
## Script: ScMaSigPro ############
##################################

go_enrichment <- function(gene_list, 
                          background,
                          ont = "BP",
                          pAdjustMethod = "BH", 
                          nterms = 5, 
                          sig.level = 0.05) {
  # Load
  suppressPackageStartupMessages(library(clusterProfiler))
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  suppressPackageStartupMessages(library(enrichplot))
  suppressPackageStartupMessages(library(AnnotationDbi))
  
  # Run GO enrichment
  ego <- enrichGO(
    gene = as.character(gene_list),
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    universe = as.character(background),
    ont = ont, # Biological Processes
    pAdjustMethod = pAdjustMethod, # Benjamini-Hochberg adjustment
    qvalueCutoff = sig.level, # Set threshold for q-value
    readable = TRUE
  ) # To show readable gene names
  
  # Check
  if (dim(ego)[1] >= 1) {
    # Calculate ema
    ego_similarity <- pairwise_termsim(ego)
    
    return(list(
      ego = ego,
      dot = dotplot(ego, showCategory = nterms),
      ema = emapplot(ego_similarity, showCategory = nterms)
    )
    )
  } else {
    return(NULL)
  }
}