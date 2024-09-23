#################################################
# Name: Priyansh Srivastava #####################
# Email: spriyansh29@gmail.com ##################
# Web: https://www.metapriyansh.com/ ############
# Script Title: Calculate Scaling Factor Table ##
#################################################

## Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(Cairo))

## Set the I/O
args <- commandArgs(trailingOnly = TRUE)
inPath_spikeInDepth <- args[1]
outPathImg <- args[2]
outPathTab <- args[3]

## load the log files of samtools
frag.list <- lapply(list.files(path = inPath_spikeInDepth, pattern = "*.seqDepth"),
  FUN = function(file) {
    return(file)
  }
)
names(frag.list) <- str_remove(frag.list, pattern = ".seqDepth")


scaleFactor <- c()
multiplier <- 10000
for (hist in names(frag.list)) {
  spikeDepth <- read.table(paste(inPath_spikeInDepth, frag.list[[hist]], sep = "/"), header = FALSE, fill = TRUE)$V1[1]

  histInfo <- str_split_1(hist, "_")[2]
  replicate <- str_split_1(hist, "_")[1]
  scaleFactor <- data.frame(scaleFactor = multiplier / spikeDepth, Histone = histInfo, Replicate = replicate) %>% rbind(scaleFactor, .)
}
scaleFactor$Histone <- as.factor(scaleFactor$Histone)

# Load Alignmnet summary
alignSummary <- read.table(paste(outPathTab, "01_Alignment_Summary.txt", sep = ""),
  header = T
)

alignSummaryNormal <- left_join(alignSummary, scaleFactor, by = c("Histone", "Replicate"))

# Calc depth
alignSummaryNormal$normDepth <- alignSummaryNormal$MappedFragNum_hg38 * alignSummaryNormal$scaleFactor

## Generate sequencing depth boxplot
fig6A <- scaleFactor %>% ggplot(aes(x = Histone, y = scaleFactor, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 20) +
  ylab("Spike-in Scalling Factor") +
  xlab("")


fig6B <- alignSummaryNormal %>% ggplot(aes(x = Histone, y = normDepth, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 20) +
  ylab("Normalization Fragment Count") +
  xlab("") +
  coord_cartesian(ylim = c(1000000, 130000000))


comb.plot <- ggarrange(fig6A, fig6B, ncol = 2, common.legend = TRUE, legend = "bottom")

# Save
ggsave(
  plot = fig6A,
  filename = "04_A_SpikeIn_Scaling_Factor.png",
  path = outPathImg,
  dpi = 800, width = 8, height = 6
)

# Save
ggsave(
  plot = fig6B,
  filename = "04_B_Normalization_Fragment_Count.png",
  path = outPathImg,
  dpi = 800, width = 8, height = 6
)

# Save
ggsave(
  plot = comb.plot,
  filename = "04_Combined_Scaling_Normalization.png",
  path = outPathImg,
  dpi = 800, width = 12, height = 8
)

# Write tables
write.table(alignSummaryNormal,
  file = paste(outPathTab, "04_Alignment_Summary_Scaling.txt", sep = ""),
  quote = F, sep = "\t", row.names = F, col.names = T
)

alignSummaryNormal$fileName <- paste0(alignSummaryNormal$Replicate, "_", alignSummaryNormal$Histone, ".bed")

# Write tables
write.table(alignSummaryNormal[, c("fileName", "scaleFactor"), drop = F],
  file = paste(outPathTab, "04_Scaling_factors.txt", sep = ""),
  quote = F, sep = "\t", row.names = F, col.names = T
)
