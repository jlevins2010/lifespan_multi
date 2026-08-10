library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(ggrepel)

# ------------------------------------------------------------------------------
# 1. Define Function
# ------------------------------------------------------------------------------
plot_gene_peak_scatter <- function(seurat_obj,
                                   gene,
                                   peak,
                                   celltype_col = "cellType_age",
                                   output_dir = NULL,
                                   save_pdf = TRUE,
                                   width = 6,
                                   height = 5.5) {
  
  if (!gene %in% rownames(seurat_obj[["RNA"]])) {
    warning(paste("Gene", gene, "not found in RNA assay."))
    return(NULL)
  }
  if (!peak %in% rownames(seurat_obj[["ATAC"]])) {
    warning(paste("Peak", peak, "not found in ATAC assay."))
    return(NULL)
  }
  
  # Set Identity
  Idents(seurat_obj) <- celltype_col
  
  # 1. Extract & Normalize RNA
  DefaultAssay(seurat_obj) <- "RNA"
  rna_df <- FetchData(seurat_obj, vars = c(gene, celltype_col)) %>%
    rename(Value = 1, CellType = 2) %>%
    group_by(CellType) %>%
    summarize(MeanRNA = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(RelRNA = MeanRNA / max(MeanRNA, na.rm = TRUE))
  
  # 2. Extract & Normalize ATAC
  DefaultAssay(seurat_obj) <- "ATAC"
  atac_df <- FetchData(seurat_obj, vars = c(peak, celltype_col)) %>%
    rename(Value = 1, CellType = 2) %>%
    group_by(CellType) %>%
    summarize(MeanATAC = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(RelATAC = MeanATAC / max(MeanATAC, na.rm = TRUE))
  
  # 3. Merge Datasets
  scatter_df <- inner_join(rna_df, atac_df, by = "CellType")
  
  # 4. Generate Scatter Plot
  p <- ggplot(scatter_df, aes(x = RelATAC, y = RelRNA, label = CellType)) +
    geom_point(aes(color = CellType), size = 3.5, show.legend = FALSE) +
    geom_text_repel(
      size = 3.5,
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.3
    ) +
    scale_x_continuous(labels = scales::percent, limits = c(-0.02, 1.05)) +
    scale_y_continuous(labels = scales::percent, limits = c(-0.02, 1.05)) +
    labs(
      title = paste("Expression vs. Accessibility:", gene),
      subtitle = paste("Peak:", peak),
      x = "Relative Accessibility (% of Max Cell Type)",
      y = "Relative Expression (% of Max Cell Type)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
    )
  
  # 5. Export as Vector PDF
  if (save_pdf && !is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    clean_peak <- gsub("[:\\-]", "_", peak)
    pdf_filename <- file.path(output_dir, paste0(gene, "_", clean_peak, "_scatter.pdf"))
    
    ggsave(
      filename = pdf_filename,
      plot = p,
      device = cairo_pdf,
      width = width,
      height = height,
      units = "in"
    )
    
    message("Saved PDF to: ", pdf_filename)
  }
  
  return(p)
}

output_dir <- "/home/levinsj/MultiOme/Analysis/Rscripts/Rplots"

#----------------- EYA1 plot
seurat_obj <- readRDS("/home/levinsj/MultiOme/MergedObjects/R_objects/seurat_outputs/paired_GWAS_seurat.rds")
target_types <- c("Fetal_PT", "Adult_PT", "Fetal_NPC", "Fetal_Int", "Pediatric_PT", "Adult_iPT")
seurat_obj$cellType_age <- paste0(seurat_obj$type, "_", seurat_obj$cellType)
seurat_obj <- subset(seurat_obj, subset = cellType_age %in% target_types)

plot_gene_peak_scatter(seurat_obj, gene = "EYA1", peak = "chr8-71737715-71738216", celltype_col = "cellType_age", output_dir = output_dir)

#------------------- PLPP3 plot
seurat_obj <- readRDS("/home/levinsj/MultiOme/MergedObjects/R_objects/seurat_outputs/paired_GWAS_seurat.rds")
target_types <- c("Fetal_Endothelium", "Adult_Endothelium", "Pediatric_Endothelium", "Fetal_NPC", "Adult_PT", "Adult_iPT")
seurat_obj$cellType_age <- paste0(seurat_obj$type, "_", seurat_obj$cellType)
seurat_obj <- subset(seurat_obj, subset = cellType_age %in% target_types)

plot_gene_peak_scatter(seurat_obj, gene = "PLPP3", peak = "chr1-56456662-56457163", celltype_col = "cellType_age", output_dir = output_dir)
	
#--------------------- SLC36A2 plot
seurat_obj <- readRDS("/home/levinsj/MultiOme/MergedObjects/R_objects/seurat_outputs/paired_GWAS_seurat.rds")
target_types <- c("Fetal_PT", "Adult_PT", "Fetal_NPC", "Fetal_Int", "Pediatric_PT", "Adult_iPT")
seurat_obj$cellType_age <- paste0(seurat_obj$type, "_", seurat_obj$cellType)
seurat_obj <- subset(seurat_obj, subset = cellType_age %in% target_types)

plot_gene_peak_scatter(seurat_obj, gene = "SLC36A2", peak = "chr5-151321550-151322051", celltype_col = "cellType_age", output_dir = output_dir)
