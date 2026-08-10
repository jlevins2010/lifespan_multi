message("loading libraries")
library(stringr)
library(Signac)
library(Seurat)
library(ggplot2)
library(openxlsx)
library(GenomicRanges)
library(reticulate)
library(Matrix)
library(scMultiMap)
library(EnsDb.Hsapiens.v86)
library(patchwork)

message()
message("libraries loaded")
gene <- "ERBB4"


region_to_plot <- GRanges(seqnames = "chr2",
			ranges = IRanges(start = 212218000, end = 212219000))
region_to_plot_ucsc <- paste0(as.character(seqnames(region_to_plot)),"-", start(region_to_plot), "-", end(region_to_plot))

uniform_fs <- 10

concatenate_seurat_meta <- function(seurat_object, col1_name, col2_name, new_col_name) {

  # We extract the first column of the resulting data frame using [, 1]
  concatenated_values <- paste(
    seurat_object[[col1_name]][, 1],
    seurat_object[[col2_name]][, 1],
    sep = "_"
  )

  # Assign the new concatenated vector back to the Seurat object
  seurat_object[[new_col_name]] <- concatenated_values

  return(seurat_object)
}


message("loading multimap files")

input_seurat <- "/home/levinsj/MultiOme/MergedObjects/R_objects/conversion_files/multiOme_paired.rds"

bigwig_files <- list(
    fetal_Ref = "/home/levinsj/MultiOme/Analysis/SALSA/Plotting_ASA/bigwig/ERBB4/Merged_NPC_Allele_A.bw",
    fetal_Alt = "/home/levinsj/MultiOme/Analysis/SALSA/Plotting_ASA/bigwig/ERBB4/Merged_NPC_Allele_G.bw",
    Int_Ref = "/home/levinsj/MultiOme/Analysis/SALSA/Plotting_ASA/bigwig/ERBB4/Merged_Int_Allele_A.bw",
    Int_Alt = "/home/levinsj/MultiOme/Analysis/SALSA/Plotting_ASA/bigwig/ERBB4/Merged_Int_Allele_G.bw",
    PT_Ref = "/home/levinsj/MultiOme/Analysis/SALSA/Plotting_ASA/bigwig/ERBB4/Merged_PT_Allele_A.bw",
    PT_Alt = "/home/levinsj/MultiOme/Analysis/SALSA/Plotting_ASA/bigwig/ERBB4/Merged_PT_Allele_G.bw",
    DCT_Ref = "/home/levinsj/MultiOme/Analysis/SALSA/Plotting_ASA/bigwig/ERBB4/Merged_DCT_CNT_Allele_A.bw",
    DCT_Alt = "/home/levinsj/MultiOme/Analysis/SALSA/Plotting_ASA/bigwig/ERBB4/Merged_DCT_CNT_Allele_G.bw"
)

message("making bigwig track")
bigwig_plot <- BigwigTrack(
        region = region_to_plot,
        bigwig = bigwig_files,
        type = "coverage",
        bigwig.scale = "common",
        ymax = 5000)
message("saving track")

ggsave(filename = paste0("/home/levinsj/MultiOme/Analysis/Rscripts/Rplots/ASA_",gene,"_plot.eps"),
        device= "eps",
        plot = bigwig_plot,
        width = 10,
        height = 12,
        units = "in",
        dpi = "retina"
        )

message("done")

obj <- readRDS(input_seurat)
obj <- concatenate_seurat_meta(
    seurat_object = obj,
    col1_name = "cellType",
    col2_name = "type",
    new_col_name = "anno_type"
)
selected_groups <- c("NPC_Fetal","Int_Fetal","PT_Fetal","PT_Pediatric","PT_Adult")

expr_plot <- ExpressionPlot(
  object = obj,
  assay = "RNA",
  slot = "counts",
  features = gene,
  idents = selected_groups,
  group.by='anno_type') + theme(text = element_text(size = uniform_fs))
message("expression plot complete")

expr_output <- paste0("/home/levinsj/MultiOme/Analysis/Rscripts/Rplots/", gene, "_expression_by_age.eps")
ggsave(filename = expr_output, plot = expr_plot, width = 16, height = 8, device = "eps")

selected_groups <- c("LOH_Adult", "DCT_CNT_Adult", "PT_Adult")

expr_plot <- ExpressionPlot(
  object = obj,
  assay = "RNA",
  slot = "counts",
  features = gene,
  idents = selected_groups,
  group.by='anno_type') + theme(text = element_text(size = uniform_fs))
message("expression plot complete")

expr_output <- paste0("/home/levinsj/MultiOme/Analysis/Rscripts/Rplots/", gene, "_expression_by_age_adult.eps")
ggsave(filename = expr_output, plot = expr_plot, width = 16, height = 8, device = "eps")

