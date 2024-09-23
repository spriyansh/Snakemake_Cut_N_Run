runDESeq2 <- function(count_table, metadata, formula, altHypothesis = "greaterAbs",
                      independentFiltering = FALSE, img_out_path, tab_out_path,
                      return_table = FALSE) {
  tryCatch(
    expr = {
      # Create Object
      dds_ob <- DESeqDataSetFromMatrix(
        countData = as.data.frame(count_table),
        colData = DataFrame(metadata),
        design = as.formula(formula)
      )

      # Estimate
      dds_ob <- DESeq(dds_ob)

      # Get results
      dds_ob <- results(dds_ob,
        independentFiltering = independentFiltering,
        altHypothesis = altHypothesis
      )

      # Create file_name
      if (grepl(x = formula, pattern = "\\+")) {
        img_name <- str_replace(str_remove(formula, pattern = "~"), pattern = "\\+", replacement = "_and_")
      } else if (grepl(x = formula, pattern = "\\:")) {
        img_name <- str_replace(str_remove(formula, pattern = "~"), pattern = "\\:", replacement = "_interact_")
      } else if (grepl(x = formula, pattern = "\\*")) {
        img_name <- str_replace(str_remove(formula, pattern = "~"), pattern = "\\*", replacement = "_and_interact_")
      } else {
        img_name <- str_remove(formula, pattern = "~")
      }

      # Create Dir
      dir.create(paste(img_out_path, "06_DESeq2_MA_Plots", sep = "/"), showWarnings = FALSE, recursive = TRUE)

      output_file_name <- paste(img_out_path, "06_DESeq2_MA_Plots", paste0(img_name, ".png"), sep = "/")
      # Open a PNG device
      #png(filename = output_file_name, width = 8, height = 6, units = "in", res = 600)
      Cairo::CairoPNG(filename = output_file_name, width = 8, height = 6, units = "in", res = 600)
      # Create MA Plot
      plotMA(dds_ob,
        main = formula,
        colSig = "#F58A53",
        colNonSig = "#15918A",
        colLine = "#EE446F"
      )

      # Close the device to save the file
      dev.off()

      # Return Data

      if (return_table) {
        return(as.data.frame(dds_ob))
      } else {
        dds_df <- as.data.frame(dds_ob)

        dds_df <- dds_df %>% rownames_to_column(
          "rownames"
        )

        # Create Dir
        dir.create(paste(tab_out_path, "07_DESeq2_Tabs", sep = "/"), showWarnings = FALSE, recursive = TRUE)
        output_file_name <- paste(tab_out_path, "07_DESeq2_Tabs", paste0(img_name, ".txt"), sep = "/")

        # Write
        write.table(dds_df,
          file = output_file_name,
          sep = "\t",
          row.names = FALSE,
          quote = FALSE
        )
        dds_df_sub  <- dds_df[dds_df$padj <= 0.05,,drop=FALSE]
        output_file_name <- paste(tab_out_path, "07_DESeq2_Tabs", paste0(img_name,"_sig", ".txt"), sep = "/")
        
        # Write
        write.table(dds_df_sub,
                    file = output_file_name,
                    sep = "\t",
                    row.names = FALSE,
                    quote = FALSE
        )
      }
    },
    error = function(e) {
      message(paste("Cannot Construct", formula, paste0("and ", e$message, ".")))
    }
  )
}
