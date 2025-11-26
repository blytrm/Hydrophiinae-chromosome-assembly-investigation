# Load Packages & Libraries
packages_needed <- c("tidyverse", "scales", "RColorBrewer")
packages_to_install <- packages_needed[!(packages_needed %in% installed.packages()[,"Package"])]

if(length(packages_to_install)) {
  install.packages(packages_to_install)
}

library(tidyverse)
library(scales)
library(RColorBrewer)

# Read the BED file (with column specification)
df_bed <- read_tsv("/Users/billytrim/Desktop/Sanders_Lab/H.major_C2_Assessment/1_ch2_coverage/raw/chr2-per_base-coverage.bed", col_names = c("chrom", "start", "end", "coverage")) 
# inspect data
print(paste("Dataset dimensions:", nrow(df_bed), ncol(df_bed)))
df_bed
head(df_bed)

# Calculating Midpoint
df_bed <- df_bed %>%
  mutate(midpoint = floor((start + end) / 2))
head(df_bed)

# Quality Control Check
cat(sum(df_bed$coverage == 0))
cat(max(df_bed$coverage))
print(summary(df_bed$coverage))

# PLOTTING

# Coverage Plot
cov <- ggplot(df_bed, aes(x = midpoint, y = coverage)) +
  geom_step(color = "steelblue") +
  labs(x = "Midpoint", y = "Coverage") +
  theme_minimal()

# Ch2 Binned Mean Coverage Plot
window_size <- 1000000
df_binned <- df_bed %>%
  mutate(bin = floor(midpoint / window_size) * window_size) %>%
  group_by(bin) %>%
  summarize(midpoint = mean(midpoint), coverage = mean(coverage))
binmeancov <- ggplot(df_binned, aes(x = midpoint, y = coverage)) +
  geom_line(color = "steelblue") +
  labs(x = "Genomic Position (midpoint)",
       y = "Mean Coverage (per bin)") +
  theme_minimal()

# Ch2 average coverage below threshold 30 (bin = 1000000)
below30 <- ggplot(df_binned, aes(x = midpoint, y = coverage)) +
  geom_line(color = "steelblue") +
  geom_hline(yintercept = 30, linetype = "dashed", color = "red") +
  labs(x = "Genomic Position (midpoint)",
       y = "Mean Coverage (per bin)",) +
  theme_minimal()

# Ch2 average coverage below 30 highlighted (bin = 1000000)
below30h <- ggplot(df_binned, aes(x = midpoint, y = coverage)) +
  geom_line(color = "steelblue") +
geom_linerange(data = subset(df_binned, coverage < 30),
               aes(x = midpoint, ymin = 0, ymax = coverage),
               color = "red", size = 1, alpha = 0.7) +
  labs(x = "Genomic Position (midpoint)",
       y = "Mean Coverage (per bin)") +
  theme_minimal()

print(cov)
print(binmeancov)
print(below30)
print(below30h)
