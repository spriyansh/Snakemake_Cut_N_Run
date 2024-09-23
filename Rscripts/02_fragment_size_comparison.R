###############################################
# Name: Priyansh Srivastava ###################
# Email: spriyansh29@gmail.com ################
# Web: https://www.metapriyansh.com/ ##########
# Script Title: Compare Mapped Fragment Size ##
###############################################

## Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggpubr))

## Set the I/O
args <- commandArgs(trailingOnly = TRUE)
inPath_hg38 <- args[1]
outPathImg <- args[2]
outPathTab <- args[3]

## load the log files of samtools
frag.list <- lapply(list.files(path = inPath_hg38, pattern = "*.txt"),
  FUN = function(file) {
    return(file)
  }
)
names(frag.list) <- str_remove(frag.list, pattern = ".txt")

# Function to generate alignment results
generatefragLenult <- function(frag.list, inPath) {
  result <- lapply(frag.list, function(file) {
    # Histone and Replicate
    split_name <- str_split(file, pattern = "_", simplify = TRUE)
    histInfo <- str_remove(split_name[2], ".txt")
    Replicate <- split_name[1]
    sampleInfo <- paste(Replicate, histInfo, sep = "_")


    # Load Log
    fragLen <- read.table(paste(inPath, file, sep = "/"), header = FALSE) %>%
      mutate(
        fragLen = V1 %>% as.numeric(), fragCount = V2
        %>% as.numeric(), Weight = as.numeric(V2) / sum(as.numeric(V2)),
        Histone = histInfo, Replicate = Replicate, sampleInfo = sampleInfo
      )
  })

  # Combine into one dataframe
  do.call(rbind, result)
}

# Get the alignment results
fragLen.hg38 <- generatefragLenult(frag.list, inPath_hg38)

# Convert to factors
fragLen.hg38$sampleInfo <- as.factor(fragLen.hg38$sampleInfo)
fragLen.hg38$Histone <- as.factor(fragLen.hg38$Histone)


## Generate the fragment size density plot (violin plot)
fig5A <- fragLen.hg38 %>% ggplot(aes(x = sampleInfo, y = fragLen, weight = Weight, fill = Histone)) +
  geom_violin(bw = 5) +
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 20) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Fragment Length") +
  xlab("")

# Save
ggsave(
  plot = fig5A,
  filename = "02_A_Fragment_Length_Violin.png",
  path = outPathImg,
  dpi = 800, width = 12, height = 6
)

fig5B <- fragLen.hg38 %>% ggplot(aes(x = fragLen, y = fragCount, color = Histone, group = sampleInfo, linetype = Replicate)) +
  geom_line(linewidth = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))

# Save
ggsave(
  plot = fig5B,
  filename = "02_B_Fragment_Length.png",
  path = outPathImg,
  dpi = 800, width = 12, height = 6
)

# Plot together
combined.plot <- ggarrange(fig5A, fig5B, ncol = 1)

# Save
ggsave(
  plot = combined.plot,
  filename = "02_Combined_Fragment_Length.png",
  path = outPathImg,
  dpi = 800, width = 12, height = 12
)

# Write tables
write.table(fragLen.hg38,
  file = paste0(outPathTab, "02_Fragment_Length.txt"),
  quote = F, sep = "\t", row.names = F, col.names = T
)
