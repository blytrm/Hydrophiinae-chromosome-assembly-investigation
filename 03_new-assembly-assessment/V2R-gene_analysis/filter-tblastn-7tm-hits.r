library(tidyverse)

cnames <- c(
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
)

scaff_order <- paste0("scaffold_", 1:1242)

# Filters: alignment length == query length, alignment >= 200aa, e-value <= 1e-10
df_tblastn <- fs::dir_ls(
    "Hydrophiinae-chromosome-quality-assessment/results/v2r-annotation",
    glob = "*outfmt6"
) |>
    read_tsv(
        col_names = cnames,
        col_types = cols(),
        id = "file"
    ) |>
    mutate(
        file = str_remove(basename(file), ".7tm.*"),
        sseqid = factor(sseqid, levels = scaff_order)
    )

# Length distribution before filtering
df_tblastn |>
    ggplot(aes(x = length, fill = file)) +
    geom_step(
        aes(colour = file),
        stat = "bin",
        bins = 50,
        direction = "mid"
    ) +
    geom_histogram(bins = 50, alpha = 0.5, position = "identity") +
    scale_y_continuous(expand = c(0.01, 0)) +
    scale_x_continuous(expand = c(0.01, 0)) +
    facet_wrap(vars(file), ncol = 2) +
    theme_classic() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14)
    )


# Keep high quality hits. Allowing the same query to hit multiple loci in the
# genome. Simply want to keep any high quality alignments.
df_tblastn <- df_tblastn |>
    filter(qstart == 1, qend == length, length >= 200, evalue <= 1e-10) |>
    arrange(file, sseqid, sstart, desc(bitscore))

# NOTE: there are many redundant hits to putative loci in this figure. The goal
# is to get an idea of what the pid distribution looks like to inform filtering.
# The peak at 60 and 70(-ish) will have counts higher than the final collapsed
# version
df_tblastn |>
    ggplot(aes(x = pident, fill = file)) +
    geom_vline(xintercept = 60) +
    geom_step(
        aes(colour = file),
        stat = "bin",
        bins = 50,
        direction = "mid"
    ) +
    geom_histogram(bins = 50, alpha = 0.5, position = "identity") +
    scale_y_continuous(expand = c(0.01, 0)) +
    scale_x_continuous(expand = c(0.01, 0)) +
    facet_wrap(vars(file), ncol = 2) +
    theme(strip.background = element_blank())

# Writing the curated TBLASTN hits to BED files
df_tblastn |>
    filter(pident >= 60) |>
    # Arranging start/end coordinates for BED format. Start always smaller than end
    mutate(
        start = ifelse(sstart > send, send, sstart),
        end = ifelse(sstart > send, sstart, send)
    ) |>
    select(sseqid, start, end) |>
    arrange(sseqid, start, end) |>
    write_tsv(
        file = "Hydrophiinae-chromosome-quality-assessment/results/v2r-annotation/hmaj-default-primary_dual-scaff_no-post-join_lowQ50-n14-correct_rounds_3-v1.7tm.tblastn.v2r.bed",
        col_names = FALSE
    )

# Bash
# for i in *.bed; do
    # BN=${i%.*}; bedtools merge -i $i > $BN.merged.bed
# done
