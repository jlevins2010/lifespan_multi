message("running GWAS04 - Exporting All Links and Merged Top Link per Loci")

# Load required libraries
library(Seurat)
library(GenomicRanges)
library(data.table)
library(scMultiMap)
library(dplyr)
library(Matrix)
library(purrr)

# --- 1. Define File Paths ---
inputMultiOme <- "/home/levinsj/MultiOme/MergedObjects/R_objects/seurat_outputs/paired_GWAS_seurat.rds"
peaks_to_test_path <- "/home/levinsj/MultiOme/LDSC/BED_files/hg38_bed/peaks_GWASloci_overlap.bed"

# Output paths for ALL links per developmental stage
output_fetal     <- "/home/levinsj/MultiOme/Analysis/GWAS/all_links_fetal_V2.csv"
output_pediatric <- "/home/levinsj/MultiOme/Analysis/GWAS/all_links_pediatric_V2.csv"
output_adult     <- "/home/levinsj/MultiOme/Analysis/GWAS/all_links_adult_V2.csv"

# Output path for SINGLE merged summary (one row per leadSNP across stages)
output_merged_top <- "/home/levinsj/MultiOme/Analysis/GWAS/top_links_merged_all_stages.csv"

# Load Seurat object
multiOme_paired <- readRDS(inputMultiOme)
print(Assays(multiOme_paired))
print(table(multiOme_paired@meta.data$type))
print(colnames(multiOme_paired@meta.data))
message("loaded multiOme_paired.rds")

DefaultAssay(multiOme_paired) <- "ATAC"

# Get annotations and gene coordinates
annot <- Signac::Annotation(object = multiOme_paired[["ATAC"]])
gene.coords <- Signac:::CollapseToLongestTranscript(ranges = annot)

# Prepare full peak coordinates object (GRanges)
peaks_granges_full <- Signac::granges(x = multiOme_paired[["ATAC"]])
peaks_granges_full$custom_id <- paste(
  seqnames(peaks_granges_full),
  start(peaks_granges_full),
  end(peaks_granges_full),
  sep = "-"
)

meta <- multiOme_paired@meta.data

multiOme_list <- list(
  "Fetal_samples"     = subset(x = multiOme_paired, cells = rownames(meta[meta$type == "Fetal" | meta$sample_final == "HK3558", ])),
  "Pediatric_samples" = subset(x = multiOme_paired, cells = rownames(meta[meta$type == "Pediatric" & meta$sample_final != "HK3558", ])),
  "Adult_samples"     = subset(x = multiOme_paired, cells = rownames(meta[meta$type == "Adult", ]))
)

# --- 2. Load and Prepare Peaks Data ---
peaks_filtered_df <- data.table::fread(peaks_to_test_path, header = FALSE)
colnames(peaks_filtered_df) <- c("chrom", "start", "end", "gene", "leadSNP")
peaks_filtered_df$peak_name <- paste(peaks_filtered_df$chrom, peaks_filtered_df$start, peaks_filtered_df$end, sep = "-")

print(head(peaks_filtered_df))
message(paste0("total peaks to test: ", nrow(peaks_filtered_df)))

# --- 3. Prepare Genes ---
distance_threshold <- 1e6
genes <- list()
for (name in names(multiOme_list)) {
  DefaultAssay(multiOme_list[[name]]) <- "RNA"
  multiOme_list[[name]] <- FindVariableFeatures(multiOme_list[[name]], assay = "RNA", nfeatures = 5000, verbose = FALSE)
  genes[[name]] <- VariableFeatures(multiOme_list[[name]])
}
all_unique_genes <- unique(unlist(genes))
message(paste0("total genes to test: ", length(all_unique_genes)))

# --- 4. Main Loop - Accumulating Results ---

fetal_results_list     <- list()
pediatric_results_list <- list()
adult_results_list     <- list()

for (i in unique(peaks_filtered_df$leadSNP)) {
  message(paste0("Testing SNP: ", i))
  
  peaks_to_test <- peaks_filtered_df[peaks_filtered_df$leadSNP == i, ]$peak_name
  top_peaks_granges <- peaks_granges_full[peaks_granges_full$custom_id %in% peaks_to_test]
  
  if (length(top_peaks_granges) == 0) next
  
  int_gene_names <- all_unique_genes[all_unique_genes %in% gene.coords$gene_name]
  if (length(int_gene_names) == 0) next
  
  peak_distance_matrix <- scMultiMap:::DistanceToTSS(
    peaks = top_peaks_granges,
    genes = gene.coords[match(int_gene_names, gene.coords$gene_name)],
    distance = distance_threshold
  )
  
  summ <- Matrix::summary(peak_distance_matrix)
  if (is.null(summ) || nrow(summ) == 0) next
  
  df_pairs <- data.frame(
    gene = colnames(peak_distance_matrix)[summ$j],
    peak = rownames(peak_distance_matrix)[summ$i]
  )
  message(paste0("total df_pairs for loci is: ", nrow(df_pairs)))
  
  # Process each group and store full results
  for (type in c("Fetal", "Pediatric", "Adult")) {
    sample_key <- paste0(type, "_samples")
    
    res <- scMultiMap::scMultiMap(
      multiOme_list[[sample_key]], 
      df_pairs, 
      bsample = "sample_final", 
      gene_assay = "RNA", 
      peak_assay = "ATAC"
    )
    
    if (!is.null(res) && nrow(res) > 0) {
      res$leadSNP <- i  # Ensure SNP ID remains attached
      
      if (type == "Fetal")     fetal_results_list[[as.character(i)]]     <- res
      if (type == "Pediatric") pediatric_results_list[[as.character(i)]] <- res
      if (type == "Adult")     adult_results_list[[as.character(i)]]     <- res
    }
  }
}

# --- 5. Helper Function to Filter Top Link per Loci & Add Column Prefixes ---
get_top_per_locus <- function(df, prefix) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  # Identify p-value or correlation/effect size column dynamically
  p_col   <- intersect(c("pval", "pvalue", "p.value", "p_val", "p_val_adj"), colnames(df))[1]
  est_col <- intersect(c("estimate", "cor", "statistic", "score"), colnames(df))[1]
  
  if (!is.na(p_col)) {
    top_df <- df %>%
      group_by(leadSNP) %>%
      slice_min(order_by = .data[[p_col]], n = 1, with_ties = FALSE) %>%
      ungroup()
  } else if (!is.na(est_col)) {
    top_df <- df %>%
      group_by(leadSNP) %>%
      slice_max(order_by = abs(.data[[est_col]]), n = 1, with_ties = FALSE) %>%
      ungroup()
  } else {
    warning(paste("Could not identify p-value/score column for", prefix, "; selecting first row per locus."))
    top_df <- df %>%
      group_by(leadSNP) %>%
      slice(1) %>%
      ungroup()
  }
  
  # Rename columns except leadSNP with group prefix (e.g. fetal_gene, adult_pval)
  colnames(top_df) <- ifelse(
    colnames(top_df) == "leadSNP", 
    "leadSNP", 
    paste0(prefix, "_", colnames(top_df))
  )
  
  return(top_df)
}

# --- 6. Combine and Save Output Dataframes ---

message("Writing output CSV files...")

# 6a. Export original full CSVs for each group
if (length(fetal_results_list) > 0) {
  df_fetal <- dplyr::bind_rows(fetal_results_list)
  write.csv(df_fetal, output_fetal, row.names = FALSE)
} else { df_fetal <- NULL }

if (length(pediatric_results_list) > 0) {
  df_pediatric <- dplyr::bind_rows(pediatric_results_list)
  write.csv(df_pediatric, output_pediatric, row.names = FALSE)
} else { df_pediatric <- NULL }

if (length(adult_results_list) > 0) {
  df_adult <- dplyr::bind_rows(adult_results_list)
  write.csv(df_adult, output_adult, row.names = FALSE)
} else { df_adult <- NULL }

# 6b. Process top link per locus for each developmental stage
top_fetal     <- get_top_per_locus(df_fetal, "fetal")
top_pediatric <- get_top_per_locus(df_pediatric, "pediatric")
top_adult     <- get_top_per_locus(df_adult, "adult")

# 6c. Merge top links across stages into a single wide data frame by leadSNP
merged_top <- list(top_fetal, top_pediatric, top_adult) %>%
  purrr::compact() %>%
  purrr::reduce(full_join, by = "leadSNP")

# 6d. Write out the final merged dataset
write.csv(merged_top, output_merged_top, row.names = FALSE)