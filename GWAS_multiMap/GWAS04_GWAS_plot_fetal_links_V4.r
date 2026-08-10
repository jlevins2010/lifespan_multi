message("Initializing packages and defining functions...")

# 1. Load Libraries
library(Signac)
library(Seurat)
library(ggplot2)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(rtracklayer)
library(scales)
library(dplyr)

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
  
  # Get TSS coordinates with respect to gene direction (Strand)
  idx <- match(links_df$gene, gene.coords$gene_name)
  relevant_genes <- gene.coords[idx]
  
  tss_coords <- ifelse(
    as.character(strand(relevant_genes)) == "-",
    end(relevant_genes),
    start(relevant_genes)
  )
  
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
    score = -log10(links_df$padj + 1e-301)
  )
  
  return(links_gr)
}

# --- Main Plotting Function ---
generate_regional_plots <- function(
    seurat_obj,
    peaks_gr,
    link_data_list,
    gwas_region,
    total_region,
    output_dir = "/home/levinsj/MultiOme/Analysis/Rscripts/Rplots",
    max,
    uniform_fs = 10
) {
  
  # --- Safeguards for input types ---
  if (is.character(peaks_gr)) {
    message("Importing BED file to GRanges...")
    peaks_gr <- rtracklayer::import.bed(peaks_gr)
  }
  if (is.character(seurat_obj)) {
    message("Loading Seurat RDS object...")
    seurat_obj <- readRDS(seurat_obj)
  }
  if (is.character(link_data_list)) {
    message("Loading link CSV files...")
    loaded_links <- list()
    for (label in names(link_data_list)) {
      loaded_links[[label]] <- read.csv(link_data_list[[label]])
    }
    link_data_list <- loaded_links
  }
  
  message(sprintf("=== Running Region: %s ===", gwas_region))
  
  # 1. Extract GWAS peaks
  region_gr <- StringToGRanges(gwas_region)
  hits <- findOverlaps(peaks_gr, region_gr, type = 'within')
  extracted_regions <- peaks_gr[queryHits(hits)]
  
  # 2. Process Links and Create Plots
  link_plots <- list()
  
  for(label in names(link_data_list)) {
    df <- link_data_list[[label]]
    all_links_gr <- create_signac_links(df, seurat_obj, region_gr, assay = "ATAC")
    
    if (!is.null(all_links_gr)) {
      regional_links <- subset(all_links_gr, padj < 0.01)
      
      if (length(regional_links) > 0) {
        Links(seurat_obj[["ATAC"]]) <- regional_links
        
        p <- LinkPlot(
          object = seurat_obj[["ATAC"]],
          region = total_region
        ) +
          scale_color_gradient(
            low = "lightgrey",
            high = "navyblue",
            limits = c(0, max),
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
  
  # 3. Assemble Track Plots
  gene_plot <- AnnotationPlot(seurat_obj[['ATAC']], region = total_region) +
    theme(text = element_text(size = uniform_fs))
  
  gwas_track <- PeakPlot(
    object = seurat_obj,
    assay = "ATAC",
    region = total_region,
    peaks = extracted_regions
  ) +
    labs(y = "GWAS Loci") +
    theme(axis.title.y = element_text(size = uniform_fs, angle = 0, vjust = 0.5))
  
  all_tracks <- link_plots
  all_tracks$GWAS_Loc <- gwas_track
  all_tracks$Genes <- gene_plot
  
  track_stack <- CombineTracks(
    plotlist = all_tracks,
    heights = c(rep(2, length(link_plots)), 0.5, 1)
  )
  
  # 4. Save Regional Track Plot
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_path <- file.path(output_dir, sprintf("%s_regional_links_V2.eps", gwas_region))
  ggsave(filename = output_path, plot = track_stack, width = 14, height = 10, device = "eps")
  message(paste("Saved Regional Plot to:", output_path))
}

# --- Shared Resource Setup ---
message("Loading shared resources...")
obj <- readRDS("/home/levinsj/MultiOme/MergedObjects/R_objects/seurat_outputs/paired_GWAS_seurat.rds")

# Add Annotations to ATAC Assay
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
annotations <- keepStandardChromosomes(annotations, pruning.mode = "tidy")
genome(annotations) <- "hg38"
Annotation(obj[['ATAC']]) <- annotations

link_data_list <- list(
  Fetal_cells     = read.csv("/home/levinsj/MultiOme/Analysis/GWAS/all_links_fetal_V2.csv"),
  Pediatric_cells = read.csv("/home/levinsj/MultiOme/Analysis/GWAS/all_links_pediatric_V2.csv"),
  Adult_cells     = read.csv("/home/levinsj/MultiOme/Analysis/GWAS/all_links_adult_V2.csv")
)

peaks_gr <- rtracklayer::import.bed("/home/levinsj/MultiOme/Analysis/GWAS/allpeaks.bed")

generate_regional_plots(
  seurat_obj      = obj,
  peaks_gr        = peaks_gr,
  link_data_list  = link_data_list,
  gwas_region     = "chr8-71735854-71763433",
  total_region    = "chr8-71000000-73000000",
  max = 60
)

generate_regional_plots(
  seurat_obj      = obj,
  peaks_gr        = peaks_gr,
  link_data_list  = link_data_list,
  gwas_region     = "chr5-151313144-151343594",
  total_region    = "chr5-151000000-151600000",
  max = 17
)

generate_regional_plots(
  seurat_obj      = obj,
  peaks_gr        = peaks_gr,
  link_data_list  = link_data_list,
  gwas_region     = "chr1-56445626-56506928",
  total_region    = "chr1-56350000-57000000",
  max = 200
)

### WNT9b
generate_regional_plots(
  seurat_obj      = obj,
  peaks_gr        = peaks_gr,
  link_data_list  = link_data_list,
  gwas_region     = "chr17-46846080-46897717",
  total_region    = "chr17-45700000-47300000",
  max = 50
)

#### Plot all gene peak combinations tested

out_dir  <- "/home/levinsj/MultiOme/Analysis/Rscripts/Rplots"
out_path <- file.path(out_dir, "significant_gene_peak_breakdown.pdf")
fetal_sig <- link_data_list$Fetal_cells$padj < 0.05
adult_sig <- link_data_list$Adult_cells$padj < 0.05

df_counts <- data.frame(
  Category = factor(
    c("Shared in Fetal & Adult", "Only Fetal", "Only Adult", "Not Significant in Either"),
    levels = c("Shared in Fetal & Adult", "Only Fetal", "Only Adult", "Not Significant in Either")
  ),
  Count = c(
    sum(fetal_sig & adult_sig, na.rm = TRUE),
    sum(fetal_sig & !adult_sig, na.rm = TRUE),
    sum(!fetal_sig & adult_sig, na.rm = TRUE),
    sum(!fetal_sig & !adult_sig, na.rm = TRUE)
  )
)

p <- ggplot(df_counts, aes(x = Category, y = Count, fill = Category)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = scales::comma(Count)), vjust = -0.5, fontface = "bold", size = 3.8) +
  scale_fill_manual(values = c(
    "Shared in Fetal & Adult"   = "#66c2a5",
    "Only Fetal"                = "#8da0cb",
    "Only Adult"                = "#fc8d62",
    "Not Significant in Either" = "#b3b3b3"
  )) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  theme_classic() +
  labs(
    title = "All Gene-Peaks Combinations tested",
    subtitle = "Fetal vs. Adult (padj < 0.05)",
    x = NULL,
    y = "Number of CRE"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.text.x = element_text(size = 9.5, face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

ggsave(
  filename = out_path,
  plot = p,
  device = cairo_pdf,  
  width = 7,
  height = 5
)

####
out_dir  <- "/home/levinsj/MultiOme/Analysis/Rscripts/Rplots"
out_path <- file.path(out_dir, "significant_CRE_breakdown.pdf")

# Specify the column name representing the unique peak identifier in your tables
peak_col <- "peak"  # Change to "peak_id", "CRE", etc., if named differently

fetal_df <- link_data_list$Fetal_cells
adult_df <- link_data_list$Adult_cells

# 1. Identify unique peaks that have AT LEAST ONE significant gene link (padj < 0.05)
fetal_sig_peaks <- unique(fetal_df[[peak_col]][!is.na(fetal_df$padj) & fetal_df$padj < 0.05])
adult_sig_peaks <- unique(adult_df[[peak_col]][!is.na(adult_df$padj) & adult_df$padj < 0.05])

# 2. Get the master universe of all unique CREs tested across both datasets
all_peaks <- unique(c(fetal_df[[peak_col]], adult_df[[peak_col]]))

# 3. Check significance per unique peak across the entire peak set
is_fetal_sig <- all_peaks %in% fetal_sig_peaks
is_adult_sig <- all_peaks %in% adult_sig_peaks

# 4. Build category counts per peak
df_counts <- data.frame(
  Category = factor(
    c("Shared in Fetal & Adult", "Only Fetal", "Only Adult", "Not Significant in Either"),
    levels = c("Shared in Fetal & Adult", "Only Fetal", "Only Adult", "Not Significant in Either")
  ),
  Count = c(
    sum(is_fetal_sig & is_adult_sig),
    sum(is_fetal_sig & !is_adult_sig),
    sum(!is_fetal_sig & is_adult_sig),
    sum(!is_fetal_sig & !is_adult_sig)
  )
)

# 5. Plotting 
p <- ggplot(df_counts, aes(x = Category, y = Count, fill = Category)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = scales::comma(Count)), vjust = -0.5, fontface = "bold", size = 3.8) +
  scale_fill_manual(values = c(
    "Shared in Fetal & Adult"   = "#66c2a5",
    "Only Fetal"                = "#8da0cb",
    "Only Adult"                = "#fc8d62",
    "Not Significant in Either" = "#b3b3b3"
  )) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  theme_classic() +
  labs(
    title = "Significant CRE Breakdown",
    subtitle = "Fetal vs. Adult (padj < 0.05 for ≥1 gene link)",
    x = NULL,
    y = "Number of CREs"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.text.x = element_text(size = 9.5, face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

ggsave(
  filename = out_path,
  plot = p,
  device = cairo_pdf,  
  width = 7,
  height = 5
)

####### Plot CRE significant by type
out_dir  <- "/home/levinsj/MultiOme/Analysis/Rscripts/Rplots"
out_path <- file.path(out_dir, "significant_loci_by_leadSNP.pdf")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

fetal_loci <- link_data_list$Fetal_cells %>%
  group_by(leadSNP) %>%
  summarise(fetal_sig = any(!is.na(padj) & padj < 0.05), .groups = "drop")

adult_loci <- link_data_list$Adult_cells %>%
  group_by(leadSNP) %>%
  summarise(adult_sig = any(!is.na(padj) & padj < 0.05), .groups = "drop")

# 3. Join by leadSNP and categorize each unique locus
locus_counts <- full_join(fetal_loci, adult_loci, by = "leadSNP") %>%
  mutate(
    fetal_sig = coalesce(fetal_sig, FALSE),
    adult_sig = coalesce(adult_sig, FALSE),
    Category  = case_when(
      fetal_sig & adult_sig   ~ "Shared in Fetal & Adult",
      fetal_sig & !adult_sig  ~ "Only Fetal",
      !fetal_sig & adult_sig  ~ "Only Adult",
      TRUE                    ~ "Not Significant in Either"
    ),
    Category = factor(Category, levels = c(
      "Shared in Fetal & Adult",
      "Only Fetal",
      "Only Adult",
      "Not Significant in Either"
    ))
  ) %>%
  count(Category, name = "Count", .drop = FALSE)

# 4. Generate ggplot bar chart
p <- ggplot(locus_counts, aes(x = Category, y = Count, fill = Category)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = scales::comma(Count)), vjust = -0.5, fontface = "bold", size = 3.8) +
  scale_fill_manual(values = c(
    "Shared in Fetal & Adult"   = "#66c2a5",
    "Only Fetal"                = "#8da0cb",
    "Only Adult"                = "#fc8d62",
    "Not Significant in Either" = "#b3b3b3"
  )) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  theme_classic() +
  labs(
    title = "Significant Gene Peaks combinations Breakdown (by leadSNP)",
    subtitle = "Locus is significant if ≥1 CRE has padj < 0.05",
    x = NULL,
    y = "Number of Unique Loci (leadSNPs)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(size = 9.5, face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# 5. Export as vector PDF
ggsave(
  filename = out_path,
  plot = p,
  device = cairo_pdf,
  width = 7,
  height = 5
)

cat("Successfully saved vector PDF to:", out_path, "\n")
