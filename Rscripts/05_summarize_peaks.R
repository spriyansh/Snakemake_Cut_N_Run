#######################################
# Name: Priyansh Srivastava ############
# Email: spriyansh29@gmail.com #########
# Web: https://www.metapriyansh.com/ ###
# Script Title: Summarize SEACR Peaks ##
########################################

## Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(GenomicRanges))

## Set the I/O
args <- commandArgs(trailingOnly = TRUE)
peakFolder <- args[1]
outPathImg <- args[2]
outPathTab <- args[3]

## load the log files of samtools
peak.list <- lapply(list.files(path = peakFolder, pattern = "*.bed"),
  FUN = function(file) {
    return(file)
  }
)
names(peak.list) <- str_remove(peak.list, pattern = ".bed")

# Empty vectors
peakN <- c()
peakWidth <- c()

# Laod one-by one
for (hist in names(peak.list)) {
  # tmp antibody
  histInfo <- str_remove_all(str_split_1(hist, "_")[2], "\\b\\.control\\b|\\.top")
  replicate <- str_split_1(hist, "_")[1]
  peakType <- str_remove_all(hist, pattern = paste0(replicate, "_", histInfo, "."))

  # Load peaks
  print(peak.list[[hist]])
  peakInfo <- read.table(paste(peakFolder, peak.list[[hist]], sep = "/"), header = FALSE, fill = TRUE) %>% mutate(width = abs(V3 - V2))
  peakN <- data.frame(peakN = nrow(peakInfo), peakType = peakType, Histone = histInfo, Replicate = replicate) %>% rbind(peakN, .)
  peakWidth <- data.frame(width = peakInfo$width, peakType = peakType, Histone = histInfo, Replicate = replicate) %>% rbind(peakWidth, .)
}

# Summarize
peakN %>% select(Histone, Replicate, peakType, peakN)

# Num Peaks
num_peaks <- ggplot(peakN, aes(x = Histone, y = peakN, fill = peakType)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Replicate ~ .) +
  coord_flip() +
  labs(
    title = "Number of Peaks by Histone and Peak Type",
    x = "Histone",
    y = "Number of Peaks"
  ) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  theme_minimal()

# Save
ggsave(
  plot = num_peaks,
  filename = "05_A_Number_of_Peaks.png",
  path = outPathImg,
  dpi = 800, width = 6, height = 6
)



########################################################################

# Define varibales
histL <- c("H3K27ac", "H3K27me3")
repL <- c("S2D", "S4D", "SUT")
peakType <- c("control", "top")
peakOverlap <- c()


for (type in peakType) {
  for (hist in histL) {
    overlap.gr <- GRanges()
    for (rep in repL) {
      peakInfo <- read.table(paste(peakFolder, paste0(rep, "_", hist, ".", type, ".bed"), sep = "/"), header = FALSE, fill = TRUE)
      peakInfo.gr <- GRanges(peakInfo$V1, IRanges(start = peakInfo$V2, end = peakInfo$V3), strand = "*")
      if (length(overlap.gr) > 0) {
        overlap.gr <- overlap.gr[findOverlaps(overlap.gr, peakInfo.gr)@from]
      } else {
        overlap.gr <- peakInfo.gr
      }
    }
    peakOverlap <- data.frame(peakReprod = length(overlap.gr), Histone = hist, peakType = type, Replicate = repL) %>% rbind(peakOverlap, .)
  }
}

peakReprod <- left_join(peakN, peakOverlap, by = c("Histone", "peakType", "Replicate")) %>% mutate(peakReprodRate = peakReprod / peakN * 100)
peakReprod %>% select(Histone, Replicate, peakType, peakN, peakReprodNum = peakReprod, peakReprodRate)

peakRepro_plot <- ggplot(peakReprod, aes(x = Histone, y = log2(peakReprodRate), fill = peakType)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Replicate ~ .) +
  coord_flip() +
  labs(
    title = "Peak Reproduction Rate by Histone and Peak Type",
    x = "Histone",
    y = "log2(Peak Reproduction Rate)"
  ) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  theme_minimal()


ggsave(
  plot = peakRepro_plot,
  filename = "05_B_Peak_Reproduction.png",
  path = outPathImg,
  dpi = 800, width = 6, height = 6
)

ggsave(
  plot = ggarrange(num_peaks, peakRepro_plot, ncol = 2),
  filename = "05_Combined_Peak_Calling.png",
  path = outPathImg,
  dpi = 800, width = 10, height = 6
)

# Write tables
write.table(peakReprod,
  file = paste0(outPathTab, "05_SEACR_peak_summary.txt"),
  quote = F, sep = "\t", row.names = F, col.names = T
)
