# Chromosome 2 (0-60mb masked region) quality analysis
# Plots coverage, soft clipping, and mapping quality

# libraries
library(Rsamtools)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(plotly)

# file paths
bam_file <- "/Users/billytrim/Desktop/Sanders_Lab/H.major_C2_Assessment/1_ch2_coverage/raw/hma.chr2.70mb.bam"
bed_file <- "/Users/billytrim/Desktop/Sanders_Lab/H.major_C2_Assessment/1_ch2_coverage/raw/mask-regions.bed"

# function to read the BED and create GenomicRanges object
read_bed_file <- function(bed_file) {
bed_data <- read.table(bed_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                       col.names = c("chr", "start", "end"))
  # convert to GR object
  GRanges(seqnames = bed_data$chr,
          ranges = IRanges(start = bed_data$start + 1, 
                           end = bed_data$end),
          region_id = paste0(bed_data$chr, ":", bed_data$start, "-", bed_data$end))
}

# function to analyse chromosome 2 with sliding windows
analyse_chromosome2 <- function(bam_file, bed_file, window_size = 50000) {
  # read masked regions
  masked_regions <- read_bed_file(bed_file)
  
  # get c2 boundaries from masked regions
  chr2_start <- min(start(masked_regions))
  chr2_end <- max(end(masked_regions))
  
  # create sliding window across c2
  window_starts <- seq(chr2_start, chr2_end - window_size + 1, by = window_size)
  window_ends <- pmin(window_starts + window_size - 1, chr2_end)
  
  # initialise results dataframe
  results <- data.frame(
    position = window_starts + window_size/2, # centre of window
    coverage = 0,
    soft_clip_freq = 0,
    avg_mapq = 0,
    total_reads = 0,
    soft_clipped_reads = 0
  )
  
  # analyse each window
  for (i in 1:length(window_starts)) {
    if (i %% 10 == 0) cat("Processing window", i, "of", length(window_starts), "\n")
    
    # window range
    window_range <- GRanges(seqnames = "chr2",
                            ranges = IRanges(start = window_starts[i], 
                                             end = window_ends[i]))
    # scan parameters
    scan_params <- ScanBamParam(
      which = window_range,
      what = c("pos", "cigar", "mapq", "flag")
    )
    
    # scan BAM file
    tryCatch({
      bam_data <- scanBam(bam_file, param = scan_params)[[1]]
      
      if (length(bam_data$pos) > 0) {
        # Calculate coverage (number of reads)
        results$coverage[i] <- length(bam_data$pos)
        results$total_reads[i] <- length(bam_data$pos)
        
        # Calculate average mapping quality (excluding unmapped reads)
        mapped_reads <- !is.na(bam_data$mapq)
        if (sum(mapped_reads) > 0) {
          results$avg_mapq[i] <- mean(bam_data$mapq[mapped_reads], na.rm = TRUE)
        }
        
        # Analyze CIGAR strings for soft clipping
        cigars <- bam_data$cigar[!is.na(bam_data$cigar)]
        if (length(cigars) > 0) {
          # Count reads with soft clipping (S in CIGAR string)
          soft_clipped <- grepl("S", cigars)
          results$soft_clipped_reads[i] <- sum(soft_clipped)
          results$soft_clip_freq[i] <- sum(soft_clipped) / length(cigars)
        }
      }
    }, error = function(e) {
      message(paste("Error processing window", i, ":", e$message))
    })
  }
  
  return(results)
}

# Main plot
# Function to create the main plot with broken axis scales
create_chromosome_plot <- function(results) {
  # Prepare data for plotting
  plot_data <- results %>%
    select(position, coverage, soft_clip_freq, avg_mapq) %>%
    pivot_longer(cols = c(coverage, soft_clip_freq, avg_mapq),
                 names_to = "metric", values_to = "value")
  
  # Create labels for metrics
  plot_data$metric_label <- factor(plot_data$metric,
                                   levels = c("coverage", "soft_clip_freq", "avg_mapq"),
                                   labels = c("Coverage", 
                                              "Soft-clipping frequency", 
                                              "Mapping quality"))
  
  # Create colour palette
  colors <- c("Coverage" = "#2E86AB", "Soft-clipping frequency" = "#A23B72", "Mapping quality" = "#F18F01")
  
  # Create the plot with faceted scales for better visualisation
  p <- ggplot(plot_data, aes(x = position/1e6, y = value, color = metric_label)) +
    geom_line(size = 0.8, alpha = 0.8) +
    geom_point(size = 0.5, alpha = 0.6) +
    facet_wrap(~metric_label, scales = "free_y", ncol = 1, strip.position = "left") +
    scale_color_manual(values = colors) +
    labs(
      x = "Position (Mb)",
      y = "Value",
      color = "Metric"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray50"),
      strip.text = element_text(size = 11, face = "bold"),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA),
      legend.position = "none"  # Remove legend since we have facets
    ) +
    scale_x_continuous(labels = number_format(accuracy = 0.1)) +
    scale_y_continuous(labels = number_format(accuracy = 0.01))
  
  return(p)
}

# Run the analysis and display the plot
results <- analyse_chromosome2(bam_file, bed_file, window_size = 50000)
plot <- create_chromosome_plot(results)
print(plot)
