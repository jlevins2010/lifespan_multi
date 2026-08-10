library(readr)
library(GenomicRanges)
library(rtracklayer)
library(Seurat)
library(Signac)
library(tidyverse)
library(ggplot2)

# Paths
outputbed_path <- "/home/levinsj/MultiOme/LDSC/BED_files/hg38_bed/GWAS_loci_hongbo.bed"
peak_bed <- "/home/levinsj/MultiOme/LDSC/BED_files/hg38_bed/allpeaks.bed" ### in hg38
input_seurat <- "/home/levinsj/MultiOme/MergedObjects/R_objects/conversion_files/multiOme_paired.rds"
GWAS_peak_overlap <- "/home/levinsj/MultiOme/LDSC/BED_files/hg38_bed/peaks_GWASloci_overlap.bed" ### in hg38

uniform_fs <- 10

# Load Data
obj <- readRDS(input_seurat)

gwas_df <- read_tsv(outputbed_path, col_names = c("chrom", "start", "end", "rsID", "score", "gene"), show_col_types = FALSE)
gwas_df$strand <- "*"
gwas_loci <- makeGRangesFromDataFrame(gwas_df, keep.extra.columns = TRUE)

gwas_peak_df <- read_tsv(GWAS_peak_overlap, col_names = c("chrom", "start", "end", "gene", "rsID"), show_col_types = FALSE)
gwas_peak_df$strand <- "*"
gwas_peak_loci <- makeGRangesFromDataFrame(gwas_peak_df, keep.extra.columns = TRUE)


# --- FIX 1: Subset the GRanges object, not the path string ---
target_locus <- gwas_loci[162]

# Create the region string
total_region <- paste0(
  seqnames(target_locus), "-",
  start(target_locus) - 100000, "-",
  end(target_locus) + 100000
)

# --- 2. Gene Annotation Plot ---
gene_plot <- AnnotationPlot(
  object = obj[['ATAC']],
  region = total_region
) + theme(text = element_text(size = uniform_fs))

# --- 3. ATAC Peaks Plot ---
peak_plot <- PeakPlot(
  object = obj,
  assay = "ATAC",
  region = total_region
) + labs(y = "ATAC Peaks")

# --- 4. GWAS Loci Plot ---
gwas_track <- PeakPlot(
  object = obj,
  assay = "ATAC",
  region = total_region,
  peaks = gwas_loci
) +
  labs(y = "GWAS\nLoci") + # Added \n for better label spacing
  theme(axis.title.y = element_text(size = uniform_fs, angle = 0, vjust = 0.5))

# --- 4. GWAS Loci Peak Plot ---
gwas_peaks <- PeakPlot(
  object = obj,
  assay = "ATAC",
  region = total_region,
  peaks = gwas_peak_loci
) +
  labs(y = "GWAS\nLoci") + # Added \n for better label spacing
  theme(axis.title.y = element_text(size = uniform_fs, angle = 0, vjust = 0.5))

# Combine
final_plot <- CombineTracks(
  plotlist = list(gwas_track, peak_plot, gwas_peaks, gene_plot),
  heights = c(1, 1, 1, 2)
)

# Save
expr_output <- "/home/levinsj/MultiOme/Analysis/Rscripts/Rplots/GWAS_overlap_test_IRX3.eps"
ggsave(filename = expr_output, plot = final_plot, width = 16, height = 8, device = "eps")

message("Done. Plots saved.")
