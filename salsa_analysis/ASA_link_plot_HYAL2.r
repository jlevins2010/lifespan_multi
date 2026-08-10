# Starting m05 regional plot of Gene-Peak links
message("Starting m05 regional plot of Gene-Peak for ASA")

# 1. Load Libraries
library(Signac)
library(Seurat)
library(ggplot2)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(rtracklayer)
library(scales)

# --- Function to convert DF to Signac-compatible Links ---
create_signac_links <- function(links_df, seurat_obj, GWAS_granges_region, assay = "ATAC") {
  annot <- Annotation(seurat_obj[[assay]])
  gene.coords <- Signac:::CollapseToLongestTranscript(ranges = annot)
  
  # Filter for genes present in annotation
  links_df <- links_df[links_df$gene %in% gene.coords$gene_name, ]
  if (nrow(links_df) == 0) return(NULL)
  
  # Filter for peaks within the GWAS region
  peak_ranges <- StringToGRanges(links_df$peak)
  overlaps <- findOverlaps(peak_ranges, GWAS_granges_region, type = "within")
  links_df <- links_df[queryHits(overlaps), ]
  
  if (nrow(links_df) == 0) return(NULL)
  
  # Parse peak coordinates
  peak_split <- do.call(rbind, strsplit(as.character(links_df$peak), split = '-'))
  p_chr   <- peak_split[, 1]
  p_start <- as.numeric(peak_split[, 2])
  p_end   <- as.numeric(peak_split[, 3])
  
  # Get TSS coordinates
  idx <- match(links_df$gene, gene.coords$gene_name)
  tss_coords <- start(gene.coords[idx])
  
  # Create GRanges object
  links_gr <- GRanges(
    seqnames = p_chr,
    ranges = IRanges(
      start = pmin(p_start, tss_coords),
      end = pmax(p_end, tss_coords)
    ),
    gene = links_df$gene,
    peak = links_df$peak,
    padj = links_df$padj,
    score = -log10(links_df$padj + 1e-301) # Use small offset to avoid Inf
  )
  
  return(links_gr)
}

# 2. Parameters
gwas_region  <- "chr3-50160500-50161100"
gene1        <- "HYAL2"
gene2        <- "TRPA1"
input_seurat <- "/home/levinsj/MultiOme/MergedObjects/R_objects/conversion_files/multiOme_paired.rds"
peak_bed     <- "/home/levinsj/MultiOme/Analysis/GWAS/allpeaks.bed"
uniform_fs   <- 10
total_region <- "chr3-49522514-51028513"
fixed_max_score <- 40  # The capped limit you requested

# 3. Load Seurat Object and Peak Data
message("Loading Seurat object...")
obj <- readRDS(input_seurat)
region_gr <- StringToGRanges(gwas_region)
peaks_gr <- import.bed(peak_bed)
hits <- findOverlaps(peaks_gr, region_gr, type = 'within')
extracted_regions <- peaks_gr[queryHits(hits)]

# 4. Annotation Setup
message("Setting up annotations...")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
new_names <- paste0("chr", seqlevels(annotations))
names(new_names) <- seqlevels(annotations)
annotations <- renameSeqlevels(annotations, new_names)
annotations <- keepStandardChromosomes(annotations, pruning.mode = "tidy")
genome(annotations) <- "hg38"
Annotation(obj[['ATAC']]) <- annotations

# 5. Load External Link Data
message("Loading link CSV files...")
files <- list(
  Fetal_cells = "/home/levinsj/MultiOme/Analysis/GWAS/all_links_fetal.csv",
  Pediatric_cells = "/home/levinsj/MultiOme/Analysis/GWAS/all_links_pediatric.csv",
  Adult_cells  = "/home/levinsj/MultiOme/Analysis/GWAS/all_links_adult.csv"
)

link_data_list <- list()
for(label in names(files)) {
  if (file.exists(files[[label]])) {
    link_data_list[[label]] <- read.csv(files[[label]])
  } else {
    message(paste("Warning: File not found -", files[[label]]))
  }
}

# 6. Process Links and Create Plots
link_plots <- list()

for(label in names(link_data_list)) {
  message(paste("Processing Condition:", label))
  
  df <- link_data_list[[label]]
  all_links_gr <- create_signac_links(df, obj, region_gr, assay = "ATAC")
  
  if (!is.null(all_links_gr)) {
    regional_links <- all_links_gr
    
    if (length(regional_links) > 0) {
      Links(obj[["ATAC"]]) <- regional_links
      
      p <- LinkPlot(
        object = obj[["ATAC"]],
        region = total_region
      ) +
      scale_color_gradient(
        low = "lightgrey",
        high = "blue",
        limits = c(0, fixed_max_score),
        oob = scales::squish
      ) +
      labs(y = label) +
      theme(axis.title.y = element_text(size = uniform_fs, angle = 0, vjust = 0.5))
      
      link_plots[[label]] <- p
    } else {
      link_plots[[label]] <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No significant links") +
        theme_void() + labs(y = label)
    }
  }
}

# 7. Build Final Genomic Tracks
message("Assembling final plot tracks...")

gene_plot <- AnnotationPlot(obj[['ATAC']], region = total_region) +
             theme(text = element_text(size = uniform_fs))

gwas_track <- PeakPlot(
  object = obj,
  assay = "ATAC",
  region = total_region,
  peaks = extracted_regions
) +
  labs(y = "GWAS Loci") +
  theme(axis.title.y = element_text(size = uniform_fs, angle = 0, vjust = 0.5))

# Combine all tracks into a named list
all_tracks <- link_plots
all_tracks$GWAS_Loc <- gwas_track
all_tracks$Genes <- gene_plot

track_stack <- CombineTracks(
  plotlist = all_tracks,
  heights = c(rep(2, length(link_plots)), 0.5, 1)
)

# 8. Save Regional Plot
output_path <- sprintf("/home/levinsj/MultiOme/Analysis/Rscripts/Rplots/%s_ASA_links.eps", gwas_region)
ggsave(filename = output_path, plot = track_stack, width = 14, height = 10, device = "eps")

message("Done. Plots saved.")
