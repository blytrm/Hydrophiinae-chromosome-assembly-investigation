# Chromosome Synteny Visualisation using circlize
# Comparing old and new Hydrophis major chromosome-level assemblies

# Install and load required packages
if (!require("circlize")) install.packages("circlize", repos = "https://cran.rstudio.com/")
if (!require("dplyr")) install.packages("dplyr", repos = "https://cran.rstudio.com/")

library(circlize)
library(dplyr)

# Functions
merge_bed <- function(df) {
  df_unique <- df %>%
    distinct(chr, start, end) %>%
    arrange(chr, start, end)
  
  if(nrow(df_unique) == 0) {
    return(data.frame(chr = character(), start = numeric(), end = numeric(), stringsAsFactors = FALSE))
  }
  
  result <- data.frame(chr = character(), start = numeric(), end = numeric(), stringsAsFactors = FALSE)
  current_chr <- df_unique$chr[1]
  current_start <- df_unique$start[1]
  current_end <- df_unique$end[1]
  
  for(i in 2:nrow(df_unique)) {
    row <- df_unique[i, ]
    if(row$chr == current_chr && row$start <= current_end) {
      current_end <- max(current_end, row$end)
    } else {
      result <- rbind(result, data.frame(chr = current_chr, start = current_start, end = current_end, stringsAsFactors = FALSE))
      current_chr <- row$chr
      current_start <- row$start
      current_end <- row$end
    }
  }
  result <- rbind(result, data.frame(chr = current_chr, start = current_start, end = current_end, stringsAsFactors = FALSE))
  return(result)
}

check_v2r_overlap <- function(chr, start, end, v2r_data) {
  chr_v2r <- v2r_data %>% filter(chr == !!chr)
  if(nrow(chr_v2r) == 0) return(FALSE)
  return(any((start <= chr_v2r$end) & (end >= chr_v2r$start)))
}

# Read data files
ref_data <- read.csv("combined_ref.csv", stringsAsFactors = FALSE)
alignment <- read.csv("filtered.csv", stringsAsFactors = FALSE)

# Prepare chromosome data
old_chr <- ref_data %>% 
  filter(assembly == "old") %>%
  select(chromosome, size) %>%
  mutate(chr = paste0("old_", chromosome))

new_chr <- ref_data %>% 
  filter(assembly == "new") %>%
  select(chromosome, size) %>%
  mutate(chr = paste0("new_", chromosome))

all_chr <- rbind(
  data.frame(chr = old_chr$chr, start = 0, end = old_chr$size),
  data.frame(chr = new_chr$chr, start = 0, end = new_chr$size)
)
chr_sizes <- setNames(all_chr$end, all_chr$chr)

# Prepare synteny links
links <- alignment %>%
  mutate(
    chr1 = paste0("new_", reference_chromosome),
    start1 = reference_start,
    end1 = reference_end,
    chr2 = paste0("old_", query_chromosome),
    start2 = query_start,
    end2 = query_end,
    identity = percent_identity
  ) %>%
  select(chr1, start1, end1, chr2, start2, end2, identity)

# Filter links
links <- links %>% 
  filter(identity >= 90) %>%
  filter((end1 - start1) >= 1000 | (end2 - start2) >= 1000) %>%
  filter(chr1 %in% names(chr_sizes) & chr2 %in% names(chr_sizes)) %>%
  mutate(
    chr1_size = chr_sizes[chr1],
    chr2_size = chr_sizes[chr2]
  ) %>%
  filter(
    start1 >= 0 & end1 <= chr1_size & 
    start2 >= 0 & end2 <= chr2_size &
    start1 < end1 & start2 < end2
  ) %>%
  select(-chr1_size, -chr2_size)

if(nrow(links) > 200000) {
  links <- links %>% 
    arrange(desc(identity)) %>%
    slice_head(n = 200000)
}

# Read and merge V2R gene positions
v2r_new_raw <- read.table("bed-v2r/new-final.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names = c("chr", "start", "end"))
v2r_old_raw <- read.table("bed-v2r/old-final.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names = c("chr", "start", "end"))

v2r_new <- merge_bed(v2r_new_raw)
v2r_old <- merge_bed(v2r_old_raw)

# Prepare V2R data for plotting
v2r_new_plot <- v2r_new %>%
  mutate(chr = paste0("new_", chr)) %>%
  filter(chr %in% all_chr$chr) %>%
  mutate(chr_size = chr_sizes[chr]) %>%
  filter(start >= 0 & end <= chr_size & start < end) %>%
  select(chr, start, end)

v2r_old_plot <- v2r_old %>%
  mutate(chr = paste0("old_", chr)) %>%
  filter(chr %in% all_chr$chr) %>%
  mutate(chr_size = chr_sizes[chr]) %>%
  filter(start >= 0 & end <= chr_size & start < end) %>%
  select(chr, start, end)

v2r_all <- rbind(v2r_new_plot, v2r_old_plot)

# Prepare links with chromosome names and colours
links <- links %>%
  mutate(
    chr1_name = gsub("new_", "", chr1),
    chr2_name = gsub("old_", "", chr2)
  )

all_chr_names <- sort(unique(c(links$chr1_name, links$chr2_name)))
chr_pair_colours <- setNames(rainbow(length(all_chr_names)), all_chr_names)

# Check V2R overlaps
links <- links %>%
  mutate(
    link_colour_base = ifelse(chr1_name %in% names(chr_pair_colours), chr_pair_colours[chr1_name], "#808080"),
    overlaps_v2r_new = mapply(function(ch, st, en) check_v2r_overlap(ch, st, en, v2r_new_plot), chr1, start1, end1),
    overlaps_v2r_old = mapply(function(ch, st, en) check_v2r_overlap(ch, st, en, v2r_old_plot), chr2, start2, end2),
    overlaps_v2r = overlaps_v2r_new | overlaps_v2r_old
  )

# Main plot
png("synteny_plot.png", width = 3000, height = 3000, res = 300)

circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = c(rep(2, nrow(new_chr)-1), 10, rep(2, nrow(old_chr)-1), 10),
  cell.padding = c(0.01, 0, 0.01, 0)
)

circos.genomicInitialize(all_chr, plotType = NULL)

chr_colours <- c(rep("#cd5c5c", nrow(new_chr)), rep("#5f9ea0", nrow(old_chr)))

circos.track(
  ylim = c(0, 1),
  panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    circos.rect(xlim[1], 0, xlim[2], 1, col = chr_colours[CELL_META$sector.numeric.index], border = NA)
    circos.text(mean(xlim), 1.5, gsub("old_|new_|chr", "", chr), cex = 0.6, facing = "clockwise", niceFacing = TRUE, font = 2)
  },
  bg.border = NA,
  track.height = 0.06
)

if(nrow(v2r_all) > 0) {
  v2r_all_with_value <- v2r_all %>% mutate(value = 0.5)
  circos.genomicTrack(
    v2r_all_with_value,
    ylim = c(0, 1),
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, ytop = 1, ybottom = 0, col = "#FF8C00", border = "#FF8C00", ...)
    },
    track.height = 0.05,
    bg.border = NA
  )
}

region1 <- links %>% select(chr = chr1, start = start1, end = end1)
region2 <- links %>% select(chr = chr2, start = start2, end = end2)

link_colour_vec <- ifelse(
  links$overlaps_v2r,
  adjustcolor("#FF8C00", alpha.f = 0.85),
  adjustcolor(links$link_color_base, alpha.f = 0.6)
)

circos.genomicLink(region1, region2, col = link_colour_vec, border = NA)

legend("topright",
  legend = c("New Assembly", "Old Assembly", "V2R Genes", "V2R Links"),
  fill = c("#cd5c5c", "#5f9ea0", "#FF8C00", "#FF8C00"),
  border = NA, bty = "n", cex = 0.6, ncol = 2)

circos.clear()
dev.off()

# Subset plot (chr2 and chrZ)
chr_subset <- c("chr2", "chrZ")
all_chr_subset <- all_chr %>% filter(gsub("old_|new_", "", chr) %in% chr_subset)
links_subset <- links %>% filter(chr1_name %in% chr_subset & chr2_name %in% chr_subset)
v2r_subset <- v2r_all %>% filter(gsub("old_|new_", "", chr) %in% chr_subset)

if(nrow(links_subset) > 0) {
  png("synteny_plot_chr2_Z.png", width = 3000, height = 3000, res = 300)
  
  circos.clear()
  circos.par(
    start.degree = 90,
    gap.degree = c(rep(2, length(chr_subset)-1), 10, rep(2, length(chr_subset)-1), 10),
    cell.padding = c(0.01, 0, 0.01, 0)
  )
  
  circos.genomicInitialize(all_chr_subset, plotType = NULL)
  
  chr_subset_colours <- c(rep("#cd5c5c", length(chr_subset)), rep("#5f9ea0", length(chr_subset)))
  
  circos.track(
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      chr = CELL_META$sector.index
      xlim = CELL_META$xlim
      circos.rect(xlim[1], 0, xlim[2], 1, col = chr_subset_colours[CELL_META$sector.numeric.index], border = NA)
      circos.text(mean(xlim), 1.5, gsub("old_|new_|chr", "", chr), cex = 0.8, facing = "clockwise", niceFacing = TRUE, font = 2)
    },
    bg.border = NA,
    track.height = 0.06
  )
  
  if(nrow(v2r_subset) > 0) {
    v2r_subset_with_value <- v2r_subset %>% mutate(value = 0.5)
    circos.genomicTrack(
      v2r_subset_with_value,
      ylim = c(0, 1),
      panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, ytop = 1, ybottom = 0, col = "#FF8C00", border = "#FF8C00", ...)
      },
      track.height = 0.05,
      bg.border = NA
    )
  }
  
  chr_subset_names <- sort(unique(c(links_subset$chr1_name, links_subset$chr2_name)))
  subset_palette <- c("paleturquoise", "red")
  chr_subset_pair_colours <- setNames(subset_palette[1:length(chr_subset_names)], chr_subset_names)
  
  links_subset <- links_subset %>%
    mutate(
      link_colour_base_subset = ifelse(chr1_name %in% names(chr_subset_pair_colours), chr_subset_pair_colours[chr1_name], "#808080"),
      overlaps_v2r_new_subset = mapply(function(ch, st, en) check_v2r_overlap(ch, st, en, v2r_new_plot), chr1, start1, end1),
      overlaps_v2r_old_subset = mapply(function(ch, st, en) check_v2r_overlap(ch, st, en, v2r_old_plot), chr2, start2, end2),
      overlaps_v2r_subset = overlaps_v2r_new_subset | overlaps_v2r_old_subset
    )
  
  region1_subset <- links_subset %>% select(chr = chr1, start = start1, end = end1)
  region2_subset <- links_subset %>% select(chr = chr2, start = start2, end = end2)
  
  link_colour_vec_subset <- ifelse(
    links_subset$overlaps_v2r_subset,
    adjustcolor("#FF8C00", alpha.f = 0.85),
    adjustcolor(links_subset$link_color_base_subset, alpha.f = 0.8)
  )
  
  circos.genomicLink(region1_subset, region2_subset, col = link_colour_vec_subset, border = NA)
  
  legend("topright",
    legend = c("New Assembly", "Old Assembly", "V2R Genes", "V2R Links"),
    fill = c("#cd5c5c", "#5f9ea0", "#FF8C00", "#FF8C00"),
    border = NA, bty = "n", cex = 0.8)
  
  circos.clear()
  dev.off()
}

cat("Plots saved: synteny_plot.png, synteny_plot_chr2_Z.png\n")
