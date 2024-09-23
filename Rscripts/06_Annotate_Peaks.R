##########################################
# Name: Priyansh Srivastava ##############
# Email: spriyansh29@gmail.com ###########
# Web: https://www.metapriyansh.com/ #####
# Script Title: Annotate Peaks ###########
##########################################

## Load libraries
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(chromVAR))
suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v75))
suppressPackageStartupMessages(library(org.Hs.eg.db))

## Set the I/O
args <- commandArgs(trailingOnly = TRUE)
control_peakFolder <- args[1]
bamFolder <- args[2]
tableOutput <- args[3]
gtfFile <- args[4]

# Create TxDb object
invisible(txdb <- makeTxDbFromGFF(file = gtfFile, format = "gtf"))

# Get Master Peak
peak.list <- lapply(list.files(path = control_peakFolder, pattern = "*.control.bed"),
  FUN = function(file) {
    return(file)
  }
)
names(peak.list) <- str_remove(peak.list, pattern = ".control.bed")

## Create Granges
mPeak <- GRanges()

## overlap with bam file to get count
for (hist in names(peak.list)) {
  peakRes <- read.table(paste(control_peakFolder, peak.list[[hist]], sep = "/"), header = FALSE, fill = TRUE)
  mPeak <- GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
}

# Master Peak

invisible(masterPeak <- GenomicRanges::reduce(mPeak))

# Annotated Peaks
invisible(masterPeak_anno <- annotatePeak(masterPeak, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db"))

# Extract
masterPeak <- masterPeak_anno@anno

# Convert to data.frame
masterPeak_df <- as.data.frame(masterPeak)

# Make counts
countMat <- matrix(NA, length(masterPeak), 2 * 3)

## load the log files of samtools
bam.list <- list.files(path = bamFolder, pattern = "*.bam", full.names = TRUE)
names(bam.list) <- str_remove(list.files(path = bamFolder, pattern = "*.bam", full.names = FALSE), ".bam")

# Remove IgG
bam.list <- bam.list[-grep("IgG", names(bam.list))]

# Create Counts
i <- 1

invisible(for (n in names(bam.list)) {
  bamFile <- bam.list[[n]]
  fragment_counts <- chromVAR::getCounts(bamFile, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
  countMat[, i] <- counts(fragment_counts)[, 1]
  i <- i + 1
})
colnames(countMat) <- str_remove(names(bam.list), pattern = ".bam")

# Create Row names
masterPeak_df$rownames <- row_vector <- paste(as.vector(masterPeak_df$SYMBOL),
  as.vector(masterPeak_df$seqnames),
  as.vector(masterPeak_df$start),
  as.vector(masterPeak_df$end),
  sep = "_"
)

# Reansfer rownames
rownames(countMat) <- masterPeak_df$rownames

# Write object
write.table(
  x = countMat,
  file = paste(tableOutput, "06_A_Count_Table.txt", sep = "/"),
  quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA
)

# Write object
write.table(
  x = masterPeak_df,
  file = paste(tableOutput, "06_B_Peak_Metadata_Table.txt", sep = "/"),
  quote = FALSE, sep = "\t", row.names = FALSE
)
