message("loading libraries")
library(Signac)
library(Seurat)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)

# --- Paths ---
# Target specific gene folder
gene_dir     <- "/home/levinsj/MultiOme/Analysis/SALSA/Plotting_ASA/bigwig/ERBB4_control/"
input_seurat <- "/home/levinsj/MultiOme/MergedObjects/R_objects/conversion_files/multiOme_paired.rds"
output_dir   <- "/home/levinsj/MultiOme/Analysis/Rscripts/Rplots/"

# --- Parameters ---
# Define the region for HYAL2 (Example: chr3:50330000-50350000)
# Update these coordinates to your specific area of interest
window <- 500
chrom <- "chr2"
start_pos <- 212237901 - window
end_pos   <- 212237901 + window

bigwig_files <- list(
    NPC_ref = "/home/levinsj/MultiOme/Analysis/SALSA/Plotting_ASA/bigwig/ERBB4_control/Merged_NPC_Allele_C.bw",
    NPC_alt = "/home/levinsj/MultiOme/Analysis/SALSA/Plotting_ASA/bigwig/ERBB4_control/Merged_NPC_Allele_T.bw",
    Int_ref = "/home/levinsj/MultiOme/Analysis/SALSA/Plotting_ASA/bigwig/ERBB4_control/Merged_Int_Allele_C.bw",
    Int_alt = "/home/levinsj/MultiOme/Analysis/SALSA/Plotting_ASA/bigwig/ERBB4_control/Merged_Int_Allele_T.bw"
)
    
current_region <- GRanges(
    seqnames = chrom,
    ranges = IRanges(start = start_pos, end = end_pos)
)

# --- Plot Tracks ---
tryCatch({
    message("Generating BigWigPlot...")
    
    bw_plot <- BigwigTrack(
        region = current_region,
        bigwig = bigwig_files,
        bigwig.scale = "common" # Ensures all tracks use the same Y-axis scale
    )

    # Save output
    bw_output <- paste0("/home/levinsj/MultiOme/Analysis/Rscripts/Rplots/ERBB4_control_ASA_bigwig.eps")
    ggsave(filename = bw_output, plot = bw_plot, width = 8, height = 12, device = "eps")

    
    message(paste("Successfully saved to:", bw_output))
}, error = function(e) {
    message("Error during plotting: ", e$message)
})
