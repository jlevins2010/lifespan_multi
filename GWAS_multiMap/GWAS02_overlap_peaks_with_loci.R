message("running GWAS02")
message("this converts ")

library(readr)
library(GenomicRanges)
library(rtracklayer)

GWAS_data <- "/home/levinsj/MultiOme/Analysis/GWAS/GWAS_loci_hongbo_hg19.csv" ### in hg19 format
outputbed <- "/home/levinsj/MultiOme/LDSC/BED_files/hg38_bed/GWAS_loci_hongbo.bed" ### in hg38 format
peak_bed <- "/home/levinsj/MultiOme/LDSC/BED_files/hg38_bed/allpeaks.bed" ### in hg38
GWAS_peak_overlap <- "/home/levinsj/MultiOme/LDSC/BED_files/hg38_bed/peaks_GWASloci_overlap.bed" ### in hg38
chain_path <- "/home/levinsj/MultiOme/Analysis/GWAS/hg19ToHg38.over.chain"

data <- read_csv(GWAS_data, col_names = FALSE)
message("data loaded")

message("parsing coordinates")
parsed_locations <- strsplit(data$X1, "[:-]")
chromosomes <- sapply(parsed_locations, function(x) x[1])
start_positions_1based <- sapply(parsed_locations, function(x) as.numeric(x[2]))
end_positions_1based <- sapply(parsed_locations, function(x) as.numeric(x[3]))

bed_data_hg19 <- data.frame(
  chrom = chromosomes,
  chromStart = start_positions_1based, # Using 1-based here for GRanges creation
  chromEnd = end_positions_1based,
  name = data$X2,
  score = data$X4,
  gene = data$X7,
  stringsAsFactors = FALSE
)

gr_hg19 <- makeGRangesFromDataFrame(bed_data_hg19,
                                    keep.extra.columns = TRUE,
                                    ignore.strand = TRUE,
                                    seqnames.field = "chrom",
                                    start.field = "chromStart",
                                    end.field = "chromEnd")

message("lifting over")
chain <- import.chain(chain_path)
curated_list <- liftOver(gr_hg19, chain)
gr_hg38 <- unlist(curated_list)
bed_data_hg38 <- as.data.frame(gr_hg38)

message("writing lifted over to bed file")
final_bed_output <- data.frame(
  chrom = bed_data_hg38$seqnames,
  chromStart = bed_data_hg38$start - 1,
  chromEnd = bed_data_hg38$end,
  name = bed_data_hg38$name,
  score = bed_data_hg38$score,
  gene = bed_data_hg38$gene
)


# Write the hg38 BED file
write.table(
  final_bed_output,
  outputbed,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

message(paste("Success! Lifted over to hg38 and wrote to:", outputbed))
message(paste("Original count:", nrow(bed_data_hg19)))
message(paste("Successfully lifted count:", nrow(final_bed_output)))

message("overlapping lifted over file within peaks file")
peaks_gr <- import.bed(peak_bed)
hits <- findOverlaps(peaks_gr, gr_hg38, type = 'within')

message("convering to bed output")
q_idx <- queryHits(hits)
s_idx <- subjectHits(hits)

overlap_bed_output <- data.frame(
  chrom = seqnames(peaks_gr)[q_idx],
  chromStart = start(peaks_gr)[q_idx],
  chromEnd = end(peaks_gr)[q_idx],
  gene = mcols(gr_hg38)$gene[s_idx],
  gwas_name = mcols(gr_hg38)$name[s_idx]
)

message(paste("Total Peaks:", length(peaks_gr)))
message(paste("Peaks in GWAS loci:", nrow(overlap_bed_output)))


message("writing output bed")

write.table(
  overlap_bed_output,
  GWAS_peak_overlap,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

message(paste("Overlap complete! Found", nrow(overlap_bed_output), "peaks in GWAS loci."))
