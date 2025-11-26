# Simple R script to plot alignment width distribution
library(ggplot2)
library(dplyr)

# File 1
file2 <- "/Users/billytrim/Desktop/Sanders_Lab/H.major_C2_Assessment/2_Hmajch2-chz/Chromosome-Chromosome_Alignments/qc_paf/Interm-raw-data/aw-maj2_HcurZ.bed"
data2 <- read.table(file2, sep = "\t", col.names = c("chr", "start", "end", "width"))
data2$sample <- "maj2_HcurZ"
data2$width_abs <- abs(data2$width)
cat("File 2 loaded:", nrow(data2), "rows\n")

# File 2
file3 <- "/Users/billytrim/Desktop/Sanders_Lab/H.major_C2_Assessment/2_Hmajch2-chz/Chromosome-Chromosome_Alignments/qc_paf/Interm-raw-data/aw-maj2_HcyaZ.bed"
data3 <- read.table(file3, sep = "\t", col.names = c("chr", "start", "end", "width"))
data3$sample <- "maj2_HcyaZ"
data3$width_abs <- abs(data3$width)
cat("File 3 loaded:", nrow(data3), "rows\n")

# File 3
file4 <- "/Users/billytrim/Desktop/Sanders_Lab/H.major_C2_Assessment/2_Hmajch2-chz/Chromosome-Chromosome_Alignments/qc_paf/Interm-raw-data/aw-maj2_NscuZ.bed"
data4 <- read.table(file4, sep = "\t", col.names = c("chr", "start", "end", "width"))
data4$sample <- "maj2_NscuZ"
data4$width_abs <- abs(data4$width)
cat("File 4 loaded:", nrow(data4), "rows\n")

# Combine all data
all_data <- rbind(data2, data3, data4)

# Check the data
cat("\nCombined data structure:\n")
str(all_data)
cat("\nColumn names:", colnames(all_data), "\n")
cat("\nSample names:", unique(all_data$sample), "\n")
cat("\nTotal rows:", nrow(all_data), "\n")

# Set colors
colors <- c("maj2_HcurZ" = "#56B4E9", 
            "maj2_HcyaZ" = "#009E73", 
            "maj2_NscuZ" = "#CC79A7")

overall_mean <- mean(all_data$width_abs)

# Create density plot with log scale and labeled vertical lines
p1 <- ggplot(all_data, aes(x = width_abs, fill = sample, color = sample)) +
  geom_density(alpha = 0.3, linewidth = 1) +
  geom_vline(xintercept = overall_mean, color = "red", linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = 5000, color = "blue", linewidth = 1, linetype = "dashed") +
  annotate("text", x = overall_mean * 0.4, y = Inf, 
           label = paste("Mean:", round(overall_mean, 0), "bp"), 
           color = "red", size = 3, vjust = 3.5) +
  annotate("text", x = 6200 * 1.2, y = Inf, 
           label = "5000 bp", 
           color = "blue", size = 3, vjust = 3.5) +
  scale_x_log10(breaks = c(100, 500, 1000, 5000, 10000, 50000, 100000),
                labels = c("100", "500", "1K", "5K", "10K", "50K", "100K")) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(x = "Alignment Width (bp)",
       y = "Density") +
  theme_minimal()

# Calculate summary statistics
summary_stats <- all_data %>%
  group_by(sample) %>%
  summarise(
    n_alignments = n(),
    mean_width = mean(width_abs),
    median_width = median(width_abs),
    min_width = min(width_abs),
    max_width = max(width_abs),
    sd_width = sd(width_abs)
  )

# Print results
cat("\nSummary Statistics:\n")
print(summary_stats)

# Display plots
print(p1)

ggsave("alignment-length-dist.png", p1)
