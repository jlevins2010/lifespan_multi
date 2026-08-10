#--- loading libraries ---#
library(Seurat)
library(GenomicRanges)
library(data.table) 
message("this script generates a bed file with all peaks from the paired Seurat Object! Running now...")

# --- User-Defined Paths ---
# Define input/output paths outside the main logic
inputMultiOme <- "/home/levinsj/MultiOme/MergedObjects/R_objects/conversion_files/multiOme_paired.rds"
peakBedFile <- "/home/levinsj/MultiOme/LDSC/BED_files/hg38_bed/allpeaks.bed"

# --------------------------------------
seurat_to_bed <- function(seurat_obj, output_path, assay_name = "ATAC") {
  atac_assay <- GetAssay(seurat_obj, assay = assay_name)
  if (!inherits(atac_assay, "ChromatinAssay")) {
    stop(paste("Assay '", assay_name, "' is not a ChromatinAssay or does not exist. Please check your assay structure."))
  }
  
  granges_peaks <- granges(atac_assay)
  
  n_peaks <- length(granges_peaks)
  if (n_peaks == 0) {
      stop("The specified assay contains 0 peaks. Cannot create a BED file.")
  }
  
  peak_names <- names(granges_peaks)
  if (is.null(peak_names) || length(peak_names) == 0) {
      peak_names <- as.character(granges_peaks)
  }
  
  bed_df <- data.frame(
    chrom = as.character(seqnames(granges_peaks)), # Chromosome
    start = start(granges_peaks),                  # Peak Start location
    end = end(granges_peaks),                      # Peak End location
    name = peak_names,                             # Peak identifier
    score = rep(0, n_peaks)                        # Placeholder
  )
  
  bed_df$start <- bed_df$start - 1
  bed_output <- bed_df[, c("chrom", "start", "end", "name", "score")]
  
  data.table::fwrite(
    bed_output,
    file = output_path,
    sep = "\t",
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE
  )
  
  cat(paste("\nSuccessfully created BED file at:", output_path, "\n"))
  cat(paste("Total peaks written:", nrow(bed_output), "\n"))
}

# --- Main Execution Block ---
if (!requireNamespace("Seurat", quietly = TRUE) ||
    !requireNamespace("GenomicRanges", quietly = TRUE) ||
    !requireNamespace("data.table", quietly = TRUE)) {
  cat("Please install the required packages: Seurat, GenomicRanges (BiocManager), and data.table.\n")
} else {
  OUTPUT_FILE <- peakBedFile
  tryCatch({
    cat(paste("Attempting to load Seurat object from:", inputMultiOme, "\n"))
    seurat_data <- readRDS(inputMultiOme)
    seurat_to_bed(seurat_data, OUTPUT_FILE, assay_name = "ATAC")
    
  }, error = function(e) {
    cat(paste("An error occurred during loading or conversion:\n"))
    cat(conditionMessage(e), "\n")
  })
}
