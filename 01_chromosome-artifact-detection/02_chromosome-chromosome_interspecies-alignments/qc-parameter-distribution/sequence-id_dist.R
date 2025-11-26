# R script to plot sequence identity distribution
library(ggplot2)
library(dplyr)
library(gridExtra)

setwd("/Users/billytrim/Desktop/Sanders_Lab/H.major_C2_Assessment/2_Hmajch2-chz/Chromosome-Chromosome_Alignments/qc_paf/Interm-raw-data")

# Reset graphics device
dev.off()
graphics.off()

# Function to safely read BED file
read_bed_safe <- function(file_path, sample_name) {
  if (file.exists(file_path) && file.size(file_path) > 0) {
    data <- read.table(file_path, sep = "\t", col.names = c("chr", "start", "end", "identity"))
    if (nrow(data) > 0) {
      data$sample <- sample_name
      cat("Loaded:", sample_name, "-", nrow(data), "alignments\n")
      return(data)
    } else {
      cat("Warning:", sample_name, "- file is empty\n")
      return(NULL)
    }
  } else {
    cat("Warning:", sample_name, "- file not found or empty\n")
    return(NULL)
  }
}

# Read files with updated species names
file1 <- "si-maj2_CgasZ.bed"
file2 <- "si-maj2_HcurZ.bed"
file3 <- "si-maj2_HcyaZ.bed"
file4 <- "si-maj2_NscuZ.bed"

# Read all files with updated species names
data1 <- read_bed_safe(file1, "H. curtus")
data2 <- read_bed_safe(file2, "H. curtus")
data3 <- read_bed_safe(file3, "H. cyan")
data4 <- read_bed_safe(file4, "N. scut")

# Combine non-empty datasets
all_data <- data.frame()
if (!is.null(data1)) all_data <- rbind(all_data, data1)
if (!is.null(data2)) all_data <- rbind(all_data, data2)
if (!is.null(data3)) all_data <- rbind(all_data, data3)
if (!is.null(data4)) all_data <- rbind(all_data, data4)

# Check if we have any data
if (nrow(all_data) == 0) {
  stop("No data found in any of the BED files!")
}

cat("\nCombined data structure:\n")
str(all_data)
cat("\nColumn names:", colnames(all_data), "\n")
cat("\nSample names:", unique(all_data$sample), "\n")
cat("\nTotal rows:", nrow(all_data), "\n")

# Set colors with updated species names
available_samples <- unique(all_data$sample)
colors <- c("H. curtus" = "#E69F00", "H. cyan" = "#56B4E9", "N. scut" = "#009E73")
colors <- colors[names(colors) %in% available_samples]

# Calculate overall mean
overall_mean <- mean(all_data$identity)

# 1. Density plot with vertical lines
p1 <- ggplot(all_data, aes(x = identity, fill = sample, color = sample)) +
  geom_density(alpha = 0.3, linewidth = 1) +
  geom_vline(xintercept = overall_mean, color = "red", linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = 80, color = "blue", linewidth = 1, linetype = "dashed") +
  annotate("text", x = overall_mean + 5, y = Inf, 
           label = paste("Mean:", round(overall_mean, 1), "%"), 
           color = "red", size = 3, vjust = 2) +
  annotate("text", x = 80 + 5, y = Inf, 
           label = "80% threshold", 
           color = "blue", size = 3, vjust = 4) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(title = "Sequence Identity Distribution Comparison",
       x = "Sequence Identity (%)",
       y = "Density") +
  theme_minimal()

# 2. Violin plot
p2 <- ggplot(all_data, aes(x = sample, y = identity, fill = sample)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.8) +
  geom_hline(yintercept = 80, color = "blue", linewidth = 1, linetype = "dashed") +
  scale_fill_manual(values = colors) +
  labs(title = "Sequence Identity Distribution by Species",
       x = "Species",
       y = "Sequence Identity (%)") +
  theme_minimal() +
  theme(legend.position = "none")

# 3. Histogram with facets
p3 <- ggplot(all_data, aes(x = identity, fill = sample)) +
  geom_histogram(bins = 30, alpha = 0.7, color = "black") +
  geom_vline(xintercept = 80, color = "blue", linewidth = 1, linetype = "dashed") +
  facet_wrap(~sample, scales = "free_y") +
  scale_fill_manual(values = colors) +
  labs(title = "Sequence Identity Histograms by Species",
       x = "Sequence Identity (%)",
       y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")

# 4. Cumulative distribution plot
p4 <- ggplot(all_data, aes(x = identity, color = sample)) +
  stat_ecdf(linewidth = 1) +
  geom_vline(xintercept = 80, color = "blue", linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = colors) +
  labs(title = "Cumulative Distribution of Sequence Identity",
       x = "Sequence Identity (%)",
       y = "Cumulative Probability") +
  theme_minimal()

# Calculate summary statistics
summary_stats <- all_data %>%
  group_by(sample) %>%
  summarise(
    n_alignments = n(),
    mean_identity = mean(identity),
    median_identity = median(identity),
    min_identity = min(identity),
    max_identity = max(identity),
    sd_identity = sd(identity),
    above_80_percent = sum(identity >= 80),
    percent_above_80 = round(sum(identity >= 80) / n() * 100, 1)
  ) %>%
  arrange(desc(mean_identity))

# Print results
cat("\nSequence Identity Summary Statistics:\n")
print(summary_stats)

# Save plots individually (this avoids the graphics state issue)
ggsave("sequence_identity_density.png", p1, width = 10, height = 6, dpi = 300)
ggsave("sequence_identity_violin.png", p2, width = 10, height = 6, dpi = 300)
ggsave("sequence_identity_histograms.png", p3, width = 12, height = 8, dpi = 300)
ggsave("sequence_identity_cumulative.png", p4, width = 10, height = 6, dpi = 300)

# Create combined plot
combined_plot <- grid.arrange(p1, p2, ncol = 2)
ggsave("sequence_identity_combined.png", combined_plot, width = 16, height = 6, dpi = 300)

# Save summary statistics to file
write.csv(summary_stats, "sequence_identity_summary_stats.csv", row.names = FALSE)
cat("\nSummary statistics saved to: sequence_identity_summary_stats.csv\n")

cat("\nAll plots saved successfully!\n")
