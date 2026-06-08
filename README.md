# Sumai 3 Introgression Detection from PAF Alignment

This repository contains scripts used to identify putative Sumai 3-derived introgression regions in wheat genomes using PAF output from minimap2.

## Description
This workflow filters whole-genome alignments between reference-s and wheat cultivars such as query-r and query-g. High-identity alignments are extracted, filtered, merged into larger genomic blocks, and visualized as genome-wide introgression ideograms.
It was used as part of the comparative genomics analysis of reference-s, query-r, and query-g wheat genome assemblies.

## Requirements
- minimap2
- awk
- R (≥ 4.0)
- ggplot2
- dplyr
- readr

## Input
- PAF file generated from minimap2 alignment
- Chromosome length table
- Optional QTL coordinate table
- Optional translocation coordinate table

## Output
- Filtered BED-like alignment file
- High-identity alignment file
- Merged introgression block tables
- Genome-wide summary tables
- Introgression ideogram plot in PNG and PDF format

  
## Usage
This script is for the whole-genome alignments that were generated using Minimap2 (alignment pipeline available in this repository: https://github.com/rrubylyn-mijan/whole-genome-and-chromosome-alignment-workflow-with-minimap2). 

---
## 1. Filter PAF alignments
```bash
awk 'BEGIN{OFS="\t"}
{
    identity="NA"

    for(i=13; i<=NF; i++){
        if($i ~ /^de:f:/){
            split($i,a,":")
            identity=(1-a[3])*100
        }
    }

    if(identity != "NA" && $12 >= 60 && $11 >= 20000 && $1 == $6){
        print $6, $8, $9, identity, $11, $12
    }
}' frj_wheat_chr_1-7.paf > query_reference.filtered.alignments.bed
```


## 2. Keep high-identity alignments
```bash
awk 'BEGIN{OFS="\t"} $4 >= 99.5 && $5 >= 20000 && $6 >= 60 {
    print $1,$2,$3,$4,$5,$6
}' query_reference.filtered.alignments.bed > query-1_reference.high_identity.alignments.bed

# Introgression detection from PAF alignment
# Description: Identifies high-identity reference-like genomic regions in query-1 or query-2
# Input: BED-like filtered alignment file from minimap2 PAF
# Output: Introgression summary tables and genome-wide ideogram
```

## 3. Visualization in R
```bash
library(ggplot2)
library(dplyr)
library(readr)

# ============================================================
# INPUT / OUTPUT
# ============================================================

blocks_file <- "query-1_reference.high_identity.alignments.bed"
out_prefix <- "query-1_reference_filtered_introgression"

# ============================================================
# FILTERING SETTINGS
# ============================================================

min_identity <- 99.99
min_aligned_bp <- 1000000
min_mapq <- 60

merge_gap <- 10000000
min_final_block_size <- 10000000

expected_reference_percent <- 6.25

# ============================================================
# CHROMOSOME ORDER
# ============================================================

chr_order <- c(
  "chr1A","chr1B","chr1D",
  "chr2A","chr2B","chr2D",
  "chr3A","chr3B","chr3D",
  "chr4A","chr4B","chr4D",
  "chr5A","chr5B","chr5D",
  "chr6A","chr6B","chr6D",
  "chr7A","chr7B","chr7D"
)

# ============================================================
# QUERY-1 CHROMOSOME LENGTHS
# ============================================================

chr_lengths <- data.frame(
  Chr = chr_order,
  Chr_length = c(
    606******, 699******, 500******,
    702******, 817******, 659******,
    759******, 857******, 627******,
    759******, 696030165, 529208703,
    804******, 724******, 585******,
    627******, 741******, 504******,
    751******, 770******, 654******
  )
)

# ============================================================
# MAJOR FHB QTL REGIONS
# ============================================================

qtl_regions <- data.frame(
  QTL = c("Fhb1", "Fhb2", "Fhb4", "Fhb5"),
  Chr = c("chr3B", "chr6B", "chr4B", "chr5A"),
  Start = c(14910167, 199411293, 566618413, 202381156),
  End = c(15291335, 329850360, 608821895, 394621872)
)

# ============================================================
# ROLLAG 2A / 5A TRANSLOCATION REGIONS
# ============================================================

translocation_regions <- data.frame(
  Translocation = c("chr5A", "chr2A"),
  Chr = c("chr2A", "chr5A"),
  Start = c(30241, 251164),
  End = c(11232272, 105638351)
)

# ============================================================
# READ ALIGNMENT FILE
# Columns: Chr Start End Identity Aligned_bp MAPQ
# ============================================================

raw_blocks <- read_tsv(
  blocks_file,
  col_names = c("Chr", "Start", "End", "Identity", "Aligned_bp", "MAPQ"),
  show_col_types = FALSE
)

raw_blocks <- raw_blocks %>%
  mutate(
    Chr = as.character(Chr),
    Start = as.numeric(Start),
    End = as.numeric(End),
    Identity = as.numeric(Identity),
    Aligned_bp = as.numeric(Aligned_bp),
    MAPQ = as.numeric(MAPQ),
    BlockSize = End - Start
  ) %>%
  filter(
    Chr %in% chr_order,
    End > Start,
    !is.na(Identity),
    !is.na(Aligned_bp),
    !is.na(MAPQ)
  )

# ============================================================
# FILTER HIGH-CONFIDENCE ALIGNMENTS
# ============================================================

filtered_blocks <- raw_blocks %>%
  filter(
    Identity >= min_identity,
    Aligned_bp >= min_aligned_bp,
    MAPQ >= min_mapq
  ) %>%
  arrange(Chr, Start, End)

# ============================================================
# MERGE NEARBY BLOCKS BY CHROMOSOME
# ============================================================

merge_one_chr <- function(df, chr_name, gap_size) {
  
  df <- df %>% arrange(Start, End)
  
  if (nrow(df) == 0) {
    return(data.frame(
      Chr = character(),
      Start = numeric(),
      End = numeric(),
      Mean_identity = numeric(),
      Total_aligned_bp = numeric(),
      Alignment_count = numeric()
    ))
  }
  
  merged <- data.frame()
  
  current_start <- df$Start[1]
  current_end <- df$End[1]
  current_identity_values <- c(df$Identity[1])
  current_aligned_bp <- df$Aligned_bp[1]
  current_count <- 1
  
  if (nrow(df) > 1) {
    for (i in 2:nrow(df)) {
      
      if (df$Start[i] <= current_end + gap_size) {
        current_end <- max(current_end, df$End[i])
        current_identity_values <- c(current_identity_values, df$Identity[i])
        current_aligned_bp <- current_aligned_bp + df$Aligned_bp[i]
        current_count <- current_count + 1
      } else {
        merged <- rbind(
          merged,
          data.frame(
            Chr = chr_name,
            Start = current_start,
            End = current_end,
            Mean_identity = mean(current_identity_values),
            Total_aligned_bp = current_aligned_bp,
            Alignment_count = current_count
          )
        )
        
        current_start <- df$Start[i]
        current_end <- df$End[i]
        current_identity_values <- c(df$Identity[i])
        current_aligned_bp <- df$Aligned_bp[i]
        current_count <- 1
      }
    }
  }
  
  merged <- rbind(
    merged,
    data.frame(
      Chr = chr_name,
      Start = current_start,
      End = current_end,
      Mean_identity = mean(current_identity_values),
      Total_aligned_bp = current_aligned_bp,
      Alignment_count = current_count
    )
  )
  
  return(merged)
}

merged_blocks <- lapply(chr_order, function(chr) {
  df_chr <- filtered_blocks %>% filter(Chr == chr)
  merge_one_chr(df_chr, chr, merge_gap)
}) %>%
  bind_rows() %>%
  mutate(BlockSize = End - Start)

major_blocks <- merged_blocks %>%
  filter(BlockSize >= min_final_block_size) %>%
  arrange(Chr, Start)

# ============================================================
# SUMMARY STATISTICS
# ============================================================

total_genome_bp <- sum(chr_lengths$Chr_length)
expected_reference_bp <- total_genome_bp * expected_reference_percent / 100

observed_bp <- sum(major_blocks$BlockSize)
observed_percent <- observed_bp / total_genome_bp * 100

summary_total <- data.frame(
  Total_query-1_genome_bp = total_genome_bp,
  Expected_reference_percent = expected_reference_percent,
  Expected_reference_bp = expected_reference_bp,
  Observed_filtered_bp = observed_bp,
  Observed_filtered_percent = observed_percent,
  Number_of_raw_alignments = nrow(raw_blocks),
  Number_of_filtered_alignments = nrow(filtered_blocks),
  Number_of_merged_blocks = nrow(merged_blocks),
  Number_of_major_blocks = nrow(major_blocks),
  Min_identity_used = min_identity,
  Min_aligned_bp_used = min_aligned_bp,
  Min_MAPQ_used = min_mapq,
  Merge_gap_used = merge_gap,
  Min_final_block_size_used = min_final_block_size
)

summary_by_chr <- major_blocks %>%
  group_by(Chr) %>%
  summarise(
    Major_blocks = n(),
    Total_introgressed_bp = sum(BlockSize),
    Mean_identity = mean(Mean_identity),
    .groups = "drop"
  ) %>%
  right_join(chr_lengths, by = "Chr") %>%
  mutate(
    Major_blocks = ifelse(is.na(Major_blocks), 0, Major_blocks),
    Total_introgressed_bp = ifelse(is.na(Total_introgressed_bp), 0, Total_introgressed_bp),
    Percent_of_chromosome = Total_introgressed_bp / Chr_length * 100
  ) %>%
  arrange(match(Chr, chr_order))

# ============================================================
# SAVE TABLES
# ============================================================

write_tsv(raw_blocks, paste0(out_prefix, "_raw_blocks.tsv"))
write_tsv(filtered_blocks, paste0(out_prefix, "_filtered_blocks.tsv"))
write_tsv(merged_blocks, paste0(out_prefix, "_merged_blocks_all.tsv"))
write_tsv(major_blocks, paste0(out_prefix, "_major_blocks.tsv"))
write_tsv(summary_total, paste0(out_prefix, "_summary_total.tsv"))
write_tsv(summary_by_chr, paste0(out_prefix, "_summary_by_chr.tsv"))
write_tsv(qtl_regions, paste0(out_prefix, "_major_FHB_QTL_regions.tsv"))
write_tsv(translocation_regions, paste0(out_prefix, "_query-1_2A_5A_translocation_regions.tsv"))

print(summary_total)
print(summary_by_chr)

# ============================================================
# PREPARE PLOT DATA
# ============================================================

chr_lengths_plot <- chr_lengths %>%
  mutate(
    Chr = factor(Chr, levels = chr_order),
    StartMb = 0,
    EndMb = Chr_length / 1e6
  )

major_blocks_mb <- major_blocks %>%
  mutate(
    Chr = factor(Chr, levels = chr_order),
    StartMb = Start / 1e6,
    EndMb = End / 1e6
  )

qtl_plot <- qtl_regions %>%
  mutate(
    Chr = factor(Chr, levels = chr_order),
    StartMb = Start / 1e6,
    EndMb = End / 1e6,
    MidMb = (StartMb + EndMb) / 2
  )

translocation_plot <- translocation_regions %>%
  mutate(
    Chr = factor(Chr, levels = chr_order),
    StartMb = Start / 1e6,
    EndMb = End / 1e6,
    MidMb = (StartMb + EndMb) / 2
  )

# ============================================================
# PLOT INTROGRESSION IDEOGRAM
# ============================================================

p <- ggplot() +
  
  geom_rect(
    data = chr_lengths_plot,
    aes(
      xmin = StartMb,
      xmax = EndMb,
      ymin = as.numeric(Chr) - 0.28,
      ymax = as.numeric(Chr) + 0.28
    ),
    fill = "#D9D9D9",
    color = NA
  ) +
  
  geom_rect(
    data = major_blocks_mb,
    aes(
      xmin = StartMb,
      xmax = EndMb,
      ymin = as.numeric(Chr) - 0.28,
      ymax = as.numeric(Chr) + 0.28
    ),
    fill = "#00b7eb",
    color = NA,
    alpha = 0.9
  ) +
  
  geom_rect(
    data = qtl_plot,
    aes(
      xmin = StartMb,
      xmax = EndMb,
      ymin = as.numeric(Chr) - 0.38,
      ymax = as.numeric(Chr) + 0.38
    ),
    fill = NA,
    color = "#D89AA7",
    linewidth = 0.35
  ) +
  
  geom_text(
    data = qtl_plot,
    aes(
      x = MidMb,
      y = as.numeric(Chr) + 0.58,
      label = QTL
    ),
    size = 3,
    fontface = "bold"
  ) +
  
  geom_rect(
    data = translocation_plot,
    aes(
      xmin = StartMb,
      xmax = EndMb,
      ymin = as.numeric(Chr) - 0.35,
      ymax = as.numeric(Chr) + 0.35
    ),
    fill = NA,
    color = "#E64B35",
    linewidth = 0.35
  ) +
  
  geom_text(
    data = translocation_plot,
    aes(
      x = EndMb,
      y = as.numeric(Chr) + 0.58,
      label = Translocation
    ),
    color = "#E64B35",
    size = 3,
    fontface = "bold"
  ) +
  
  scale_y_continuous(
    breaks = seq_along(chr_order),
    labels = chr_order,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  scale_x_continuous(expand = c(0, 0)) +
  
  labs(
    title = "Putative reference-like regions in query-1 with 2A/5A translocation",
    subtitle = paste0(
      "Identity ≥ ", min_identity,
      "%; aligned length ≥ ", min_aligned_bp / 1e6,
      " Mb; merge gap ≤ ", merge_gap / 1e6,
      " Mb; final block ≥ ", min_final_block_size / 1e6,
      " Mb; observed = ", round(observed_percent, 2),
      "%, expected = ", expected_sumai3_percent, "%"
    ),
    x = "Genomic position (Mb)",
    y = "Chromosome"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 14),
    axis.line = element_line(color = "black")
  )

ggsave(
  paste0(out_prefix, "_filtered_ideogram_with_translocation.png"),
  p,
  width = 14,
  height = 8,
  dpi = 300
)

ggsave(
  paste0(out_prefix, "_filtered_ideogram_with_translocation.pdf"),
  p,
  width = 14,
  height = 8
)

print(p)
```

## Notes
The first awk script extracts same-chromosome alignments from the minimap2 PAF file.
Sequence identity is calculated from the de:f: tag in the PAF file.
The R script applies stricter filters to identify large, high-confidence Sumai 3-like genomic blocks.
The final plot shows putative introgression blocks, major FHB QTL regions, and the Rollag 2A/5A translocation.

## Example Output
query-1_reference_filtered_introgression_filtered_ideogram_with_translocation.png
query-1_reference_filtered_introgression_filtered_ideogram_with_translocation.pdf
query-1_reference_filtered_introgression_major_blocks.tsv
query-1_reference_filtered_introgression_summary_total.tsv
query-1_reference_filtered_introgression_summary_by_chr.tsv

Maintainer:

Rubylyn Mijan
