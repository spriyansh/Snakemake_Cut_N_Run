#######################################
# Name: Priyansh Srivastava ###########
# Email: spriyansh29@gmail.com ########
# Web: https://www.metapriyansh.com/ ##
# Script Title: Summarize Alignments ##
#######################################

## Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggpubr))

## Set the I/O
args <- commandArgs(trailingOnly = TRUE)
inPath_hg38 <- args[1]
inPath_spikeIn <- args[2]
outPathImg <- args[3]
outPathTab <- args[4]

## load the log files of samtools
log.list <- lapply(list.files(path = inPath_hg38, pattern = "*.log"),
  FUN = function(file) {
    return(file)
  }
)
names(log.list) <- str_remove(log.list, pattern = ".log")

# Function to generate alignment results
generateAlignResult <- function(log.list, inPath) {
  result <- lapply(log.list, function(file) {
    # Load Log
    alignRes <- read.table(paste(inPath, file, sep = "/"), header = FALSE, fill = TRUE)

    # Get alignment rate
    alignRate <- substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6])) - 1)

    # Histone and Replicate
    split_name <- str_split(file, pattern = "_", simplify = TRUE)
    histInfo <- str_remove(split_name[2], ".log")
    Replicate <- split_name[1]

    # Create df
    data.frame(
      Histone = histInfo, Replicate = Replicate,
      SequencingDepth = as.numeric(alignRes$V1[1]),
      MappedFragNum = as.numeric(alignRes$V1[4]) + as.numeric(alignRes$V1[5]),
      AlignmentRate = as.numeric(alignRate)
    )
  })

  # Combine into one dataframe
  do.call(rbind, result)
}

# Get the alignment results
alignResult.hg38 <- generateAlignResult(log.list, inPath_hg38)
alignResult.spikeIn <- generateAlignResult(log.list, inPath_spikeIn)

# Rename columns for spikein the merge
colnames(alignResult.spikeIn)[4:5] <- c("MappedFragNum_spikeIn", "AlignmentRate_spikeIn")
colnames(alignResult.hg38)[4:5] <- c("MappedFragNum_hg38", "AlignmentRate_hg38")

# Left join
alignSummary <- left_join(alignResult.hg38, alignResult.spikeIn,
  by = c("Histone", "Replicate", "SequencingDepth")
)

## Create Visual Summary
fig3A <- alignSummary %>% ggplot(aes(x = Histone, y = SequencingDepth / 1000000, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Sequencing Depth per Million") +
  xlab("") +
  ggtitle("A. Sequencing Depth")

# Save
ggsave(
  plot = fig3A,
  filename = "01_A_Sequencing_Depth.png",
  path = outPathImg,
  dpi = 800, width = 6, height = 6
)

fig3B <- alignSummary %>% ggplot(aes(x = Histone, y = MappedFragNum_hg38 / 1000000, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Mapped Fragments per Million") +
  xlab("") +
  ggtitle("B. Alignable Fragment (hg38)")

# Save
ggsave(
  plot = fig3B,
  filename = "01_B_Alignable_Fragment_hg38.png",
  path = outPathImg,
  dpi = 800, width = 6, height = 6
)

fig3C <- alignSummary %>% ggplot(aes(x = Histone, y = AlignmentRate_hg38, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% of Mapped Fragments") +
  xlab("") +
  ggtitle("C. Alignment Rate (hg38)")

# Save
ggsave(
  plot = fig3C,
  filename = "01_C_Alignment_Rate_hg38.png",
  path = outPathImg,
  dpi = 800, width = 6, height = 6
)

fig3D <- alignSummary %>% ggplot(aes(x = Histone, y = AlignmentRate_spikeIn, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Spike-in Alignment Rate") +
  xlab("") +
  ggtitle("D. Alignment Rate (E.coli)")

# Save
ggsave(
  plot = fig3D,
  filename = "01_D_Alignment_Rate_SpikeIn.png",
  path = outPathImg,
  dpi = 800, width = 6, height = 6
)

# Plot together
combined.plot <- ggarrange(fig3A, fig3B, fig3C, fig3D, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

# Save
ggsave(
  plot = combined.plot,
  filename = "01_Combined_Alignment_Summary.png",
  path = outPathImg,
  dpi = 800, width = 12, height = 12
)

# Write tables
write.table(alignSummary,
  file = paste0(outPathTab, "01_Alignment_Summary.txt"),
  quote = F, sep = "\t", row.names = F, col.names = T
)
