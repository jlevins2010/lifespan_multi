message("--- Starting Integrated scMultiMap Analysis Workflow (Human - Full Sample) ---")

library(Signac)
library(Seurat)
library(ggplot2)
library(openxlsx)
library(GenomicRanges)
library(reticulate)
library(Matrix)
library(scMultiMap)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(forcats)
library(tibble)

set.seed(42)

# --- 2. Configuration Setup ---
distance_threshold <- 5e5
chunk_size         <- 250000  # Number of pairs to process per batch to prevent RAM overflow
nGenes             <- 1000     

# Path to Seurat multiome RDS and DAR CSV file
multiomeRDS <- "/ceph/projects/ksusztak/nephrobase_1/levinsj/MultiOme/MergedObjects/R_objects/seurat_outputs/paired_PT_seurat.rds"
dar_file    <- "/vast/home/l/levinsj/DARs/Adult_PT_vs_All_Other_Adult_DARs.csv"

base_output_dir <- "/vast/home/l/levinsj/output_csv/"
dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)

merged_output          <- paste0(base_output_dir, "mergedMultiMapoutput_adultPTlinks_human.csv")
merged_output_filtered <- paste0(base_output_dir, "mergedMultiMapoutput_adultPTlinks_filtered_human.csv")

# Helper function to run scMultiMap in batches over gene-peak pairs
run_scMultiMap_batched <- function(seurat_obj, pairs_df, batch_size = chunk_size, ...) {
  n_pairs <- nrow(pairs_df)
  n_chunks <- ceiling(n_pairs / batch_size)
  results_list <- list()
  
  for (i in seq_len(n_chunks)) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, n_pairs)
    message(sprintf("    --> Processing pair batch %d/%d (%d to %d pairs)...", i, n_chunks, start_idx, end_idx))
    
    chunk_pairs <- pairs_df[start_idx:end_idx, ]
    
    res_chunk <- tryCatch({
      scMultiMap(seurat_obj, chunk_pairs, ...)
    }, error = function(e) {
      warning(paste("Error in batch", i, ":", e$message))
      return(NULL)
    })
    
    if (!is.null(res_chunk) && nrow(res_chunk) > 0) {
      results_list[[i]] <- res_chunk
    }
    
    gc() # Clean up memory after each batch
  }
  
  if (length(results_list) > 0) {
    return(bind_rows(results_list))
  } else {
    return(NULL)
  }
}

# --- 3. Load Seurat Object & Filter Peaks from DAR File ---
message("--- Loading Seurat Object & Selecting DAR-Filtered Peaks ---")
multi <- readRDS(multiomeRDS)

DefaultAssay(multi) <- "ATAC"

# Load DAR results and filter for padj < 0.01 and log2FoldChange > 1
message(paste0("Loading DAR file: ", dar_file))
dar_df <- read.csv(dar_file, row.names = 1)

dar_filtered <- dar_df %>%
  rownames_to_column(var = "peak") %>%
  filter(!is.na(padj) & padj < 0.01 & log2FoldChange > 1)

message(paste0("Number of DAR peaks passing threshold (padj < 0.01, log2FC > 1): ", nrow(dar_filtered)))

# Convert filtered peak strings (e.g., 'chr12:6215897-6216398') into GRanges using Signac
peaks_granges_full <- Signac::StringToGRanges(dar_filtered$peak, sep = c(":", "-"))

# Assign custom IDs to the GRanges object using standard hyphenated format ('chr-start-end')
peaks_granges_full$custom_id <- paste(
  seqnames(peaks_granges_full),
  start(peaks_granges_full),
  end(peaks_granges_full),
  sep = "-"
)

all_peaks <- peaks_granges_full$custom_id
message(paste0("Total peaks retained for scMultiMap: ", length(all_peaks)))

gene_list_collector <- list()

comparison_targets_genes <- list(
  list(target_condition = (multi$type == "Adult") & (multi$cellType == "PT"), name = "Adult_PT")
)

for (comp in comparison_targets_genes) {
  target_condition <- comp$target_condition
  condition_name   <- comp$name
  
  cells_to_keep <- colnames(multi)[target_condition]
  subset_cells  <- subset(x = multi, cells = cells_to_keep)
  
  rna_counts <- tryCatch(
    Seurat::GetAssayData(subset_cells, assay = "RNA", layer = "counts"),
    error = function(e) Seurat::GetAssayData(subset_cells, assay = "RNA", slot = "counts")
  )
  
  top_genes  <- names(sort(Matrix::rowSums(rna_counts), decreasing = TRUE))[1:min(nGenes, nrow(rna_counts))]
  gene_list_collector[[condition_name]] <- top_genes
}

all_unique_genes <- unique(unlist(gene_list_collector))
message(paste0("Number of unique genes: ", length(all_unique_genes)))

# --- 4. Generate Distance-Filtered Gene-Peak Pairs ---
message("--- Generating Distance-Filtered Gene-Peak Pairs ---")

annot <- Signac::Annotation(object = multi[["ATAC"]])
gene.coords <- Signac:::CollapseToLongestTranscript(ranges = annot)

int_gene_names <- all_unique_genes[all_unique_genes %in% gene.coords$gene_name]

# Calculates TSS distances only using the DAR-filtered peaks GRanges object
peak_distance_matrix <- scMultiMap:::DistanceToTSS(
  peaks = peaks_granges_full,
  genes = gene.coords[match(int_gene_names, gene.coords$gene_name)],
  distance = distance_threshold
)

summ <- Matrix::summary(peak_distance_matrix)

df_pairs <- data.frame(
  gene = colnames(peak_distance_matrix)[summ$j],
  peak = rownames(peak_distance_matrix)[summ$i]
)

message(paste0("Number of gene-peak pairs within 1MB: ", nrow(df_pairs)))

# --- 5. Configuration for scMultiMap Conditions ---
comparison_targets <- list(
  list(target_condition = (multi$type == "Fetal") & (multi$cellType == "NPC"), name = "Fetal_NPC"),
  list(target_condition = (multi$type == "Fetal") & (multi$cellType == "Int"), name = "Fetal_Int"),
  list(target_condition = (multi$type == "Fetal") & (multi$cellType %in% c("iPT", "PT")), name = "Fetal_PT"),
  list(target_condition = (multi$type == "Pediatric") & (multi$cellType == "PT"), name = "Pediatric_PT"),
  list(target_condition = (multi$type == "Adult") & (multi$cellType == "PT"), name = "Adult_PT")
)

# --- 6. Run scMultiMap Across Conditions (Batched) ---
message("--- Running scMultiMap (Batched Execution) ---")
ref_condition <- (multi$type == "Adult") & !(multi$cellType %in% c("iPT", "PT"))

for (comp in comparison_targets) {
  target_condition   <- comp$target_condition
  comp_name          <- comp$name  
  combined_condition <- target_condition | ref_condition
  
  cells_to_keep      <- colnames(multi)[combined_condition]
  subset_cells       <- subset(x = multi, cells = cells_to_keep)
  
  message(paste0("Running ", comp_name, " on ", length(Cells(subset_cells)), " cells..."))
  
  # Run scMultiMap in chunks to avoid memory limits
  res <- run_scMultiMap_batched(
    seurat_obj = subset_cells, 
    pairs_df = df_pairs, 
    batch_size = chunk_size,
    gene_assay = "RNA", 
    peak_assay = "ATAC"
  )
  
  output_file <- paste0(base_output_dir, "adultPTlinks_", comp_name, "_multiMapoutput_human.csv")
  
  if (!is.null(res) && nrow(res) > 0) {
    message(paste0("Total pairs processed for ", comp_name, ": ", nrow(res)))
    res$target_cellType <- comp_name
    write.csv(res, file = output_file, row.names = FALSE)
  } else {
    message(paste0("No results generated for ", comp_name, ". Writing empty template CSV."))
    dummy_df <- data.frame(
      gene = character(), peak = character(), pval = numeric(), 
      padj = numeric(), covar = numeric(), cor = numeric(), 
      test_stat = numeric(), target_cellType = character()
    )
    write.csv(dummy_df, file = output_file, row.names = FALSE)
  }
  
  rm(subset_cells, cells_to_keep, combined_condition, res)
  gc()
}

rm(multi)
gc()

# --- 7. Merge & Filter Results ---
desired_order <- c("Fetal_NPC", "Fetal_Int", "Fetal_PT", "Pediatric_PT", "Adult_PT")
cellTypes <- c("Adult_PT", "Pediatric_PT", "Fetal_PT", "Fetal_Int", "Fetal_NPC")
filenames <- map_chr(cellTypes, ~ str_c(base_output_dir, "adultPTlinks_", .x, "_multiMapoutput_human.csv"))

existing_files <- filenames[file.exists(filenames)]
file_list      <- set_names(map(existing_files, read.csv), existing_files)

processed_list <- imap(file_list, ~ {
  cell_type_name <- str_extract(.y, paste(cellTypes, collapse = "|"))
  
  .x %>%
    mutate(gene_peak = paste(gene, peak, sep = "_")) %>%
    rename(!!sym(paste0("padj_", cell_type_name)) := padj) %>%
    rename(!!sym(paste0("covar_", cell_type_name)) := covar) %>%
    rename(!!sym(paste0("cor_", cell_type_name)) := cor) %>%
    rename(!!sym(paste0("testStat_", cell_type_name)) := test_stat) %>%
    select(-any_of(c("pval", "gene", "peak", "target_cellType")))
})

final_merged_df <- processed_list %>% reduce(full_join, by = "gene_peak")

final_merged_df <- final_merged_df %>%
  separate(col = gene_peak, into = c("gene", "peak"), sep = "_", remove = FALSE) %>%
  select(gene_peak, gene, peak, everything())

write.csv(final_merged_df, file = merged_output, row.names = FALSE)

filtered_data <- final_merged_df %>%
  filter(if_any(.cols = starts_with("padj_"), .fns = ~ !is.na(.x) & .x < 0.05))

write.csv(filtered_data, file = merged_output_filtered, row.names = FALSE)

# --- 8. Plots ---
plot_count_data <- function(df, metric, cutoff, compare_fn, title_prefix, file_suffix, y_label) {
  metric_cols <- str_c(metric, "_", cellTypes)
  
  count_df <- df %>%
    select(any_of(metric_cols)) %>%
    summarise(across(
      .cols = everything(),
      .fns = ~ sum(compare_fn(.x, cutoff), na.rm = TRUE),
      .names = "{.col}_count"
    )) %>%
    pivot_longer(
      cols = everything(),
      names_to = "Cell_Type",
      values_to = "Count"
    ) %>%
    mutate(
      Cell_Type = gsub(paste0(metric, "_"), "", Cell_Type),
      Cell_Type = gsub("_count", "", Cell_Type),
      Cell_Type = fct_relevel(Cell_Type, desired_order)
    )
  
  sig_plot <- ggplot(count_df, aes(x = Cell_Type, y = Count, fill = Cell_Type)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = str_c(title_prefix, " (", metric, cutoff, ")"),
      x = "Cell Type",
      y = y_label
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  
  ggsave(
    filename = str_c("/vast/home/l/levinsj/plots/multiOme/adult_enhancers_human_", file_suffix, ".eps"),
    device = cairo_ps,
    plot = sig_plot,
    width = 8,
    height = 10,
    units = "in"
  )
}

plot_count_data(
  df = final_merged_df,
  metric = "padj",
  cutoff = 0.05,
  compare_fn = `<`,
  title_prefix = "Count of Significant Links",
  file_suffix = "pvalue",
  y_label = "Count of Significant Gene/Peak Links"
)

plot_count_data(
  df = final_merged_df,
  metric = "cor",
  cutoff = 0.1,
  compare_fn = `>`,
  title_prefix = "Count of Links with Positive Correlation",
  file_suffix = "pos_corr",
  y_label = "Count of Gene/Peak Links"
)

# --- 9. Correlation Conservation Scatter Plots ---
message("--- Generating Adult PT Correlation Conservation Scatter Plots ---")

# Subset for significant Adult PT links
sig_adult_df <- final_merged_df %>%
  filter(!is.na(padj_Adult_PT) & padj_Adult_PT < 0.05)

# Pivot long to pair Adult PT correlations against each developmental stage
scatter_df <- sig_adult_df %>%
  mutate(cor_Adult_PT_ref = cor_Adult_PT) %>%
  pivot_longer(
    cols = starts_with("cor_") & !matches("cor_Adult_PT_ref"),
    names_to = "Comparison_Group",
    values_to = "Correlation_Other"
  ) %>%
  rename(cor_Adult_PT = cor_Adult_PT_ref) %>%
  mutate(
    Comparison_Group = str_remove(Comparison_Group, "^cor_")
  )

# Check if scatter_df contains actual data
if (nrow(scatter_df) == 0 || all(is.na(scatter_df$cor_Adult_PT))) {
  stop("scatter_df is empty or contains only NA values. Fix the upstream matrix processing step before plotting.")
}

# 1. Reorder factor levels
scatter_df <- scatter_df %>%
  mutate(
    Comparison_Group = factor(
      Comparison_Group, 
      levels = c("Fetal_NPC", "Fetal_Int", "Fetal_PT", "Pediatric_PT", "Adult_PT")
    )
  )

# 2. Calculate R^2 safely using Pearson correlation (avoids lm() perfect fit warnings)
r2_labels <- scatter_df %>%
  filter(!is.na(cor_Adult_PT) & !is.na(Correlation_Other)) %>%
  group_by(Comparison_Group) %>%
  summarize(
    r2_label = if (n() >= 2) {
      r_val <- cor(cor_Adult_PT, Correlation_Other, use = "complete.obs")
      sprintf("R^2 == %.3f", r_val^2)
    } else {
      "R^2 == NA"
    },
    .groups = "drop"
  )

# 3. Build plot (added formula = y ~ x to silence geom_smooth message)
np_cor_scatter <- ggplot(scatter_df, aes(x = cor_Adult_PT, y = Correlation_Other)) +
  geom_point(alpha = 0.2, size = 1, color = "#2b5c8f") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "firebrick", linewidth = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 0.7, se = TRUE) +
  geom_text(
    data = r2_labels,
    aes(x = -Inf, y = Inf, label = r2_label),
    hjust = -0.2, vjust = 1.5,
    parse = TRUE,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Comparison_Group, ncol = 2) +
  labs(
    title = "Conservation of Adult PT Correlations Across Developmental Groups",
    subtitle = paste0("Filtered for significant Adult PT links (padj < 0.05, n = ", nrow(sig_adult_df), ")"),
    x = "Adult PT Correlation (r)",
    y = "Comparison Group Correlation (r)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray92"),
    axis.title = element_text(face = "bold")
  )

# 4. Save plot
ggsave(
  filename = "/vast/home/l/levinsj/plots/multiOme/adult_enhancers_human_cor_conservation.eps",
  device = cairo_ps,
  plot = np_cor_scatter,
  width = 9,
  height = 8,
  units = "in"
)
