args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript mutation_type_figures.R <data_directory> <figure_directory>")
}

data_directory   <- args[1]
figure_directory <- args[2]

# ensure figure directory exists
if (!dir.exists(data_directory)) {
  stop(sprintf("Data directory not found: '%s'", data_directory))
}

library(dplyr)
library(tidyr)
library(Biostrings)
library(GenomicRanges)
library(MutationalPatterns)
library(ggplot2)


### helper functions
add_ref_alt <- function(gr) {
  nm <- names(gr)
  if (is.null(nm)) {
    stop("`gr` must have names like 'chr1_3892395_C_A' in names(gr)")
  }
  
  parts <- strsplit(nm, "_", fixed = TRUE)
  lens  <- lengths(parts)
  
  bad   <- which(lens < 4)
  if (length(bad)) {
    stop("The following indices in names(gr) don't have ≥4 fields: ",
         paste(bad, collapse = ", "))
  }
  
  # extract the last two fields
  REF <- vapply(parts, function(x) x[length(x) - 1], character(1))
  ALT <- vapply(parts, function(x) x[length(x)],     character(1))
  
  mcols(gr)$REF <- as.character(REF)
  mcols(gr)$ALT <- as.character(ALT)
  
  return(gr)
}

plot_mutation_aggregate_percentages <- function(type_occurrences_df) {
  
  if (any(duplicated(colnames(type_occurrences_df)))) {
    stop("Column names (mutation types) must be unique.")
  }
  
  mutation_sums <- colSums(type_occurrences_df, na.rm = TRUE)
  
  mutation_df <- data.frame(
    Mutation = names(mutation_sums),
    Count = as.numeric(mutation_sums),
    stringsAsFactors = FALSE
  )
  
  ct_cpg <- mutation_df$Count[mutation_df$Mutation == "C>T at CpG"]
  ct_other <- mutation_df$Count[mutation_df$Mutation == "C>T other"]
  
  # remove C>T subtypes 
  mutation_df <- mutation_df %>%
    filter(!Mutation %in% c("C>T at CpG", "C>T other"))
  
  # new data frame for C>T with its subtypes
  ct_df <- data.frame(
    Mutation = "C>T",
    Subtype = c("C>T at CpG", "C>T other"),
    Count = c(ct_cpg, ct_other),
    stringsAsFactors = FALSE
  )
  
  total_mutations <- sum(mutation_df$Count, ct_df$Count)
  
  mutation_df <- mutation_df %>%
    mutate(Percentage = (Count / total_mutations) * 100) %>%
    select(Mutation, Percentage)
  
  ct_df <- ct_df %>%
    mutate(Percentage = (Count / total_mutations) * 100)
  
  other_df <- mutation_df %>%
    mutate(Subtype = Mutation)
  
  plot_df <- bind_rows(other_df, ct_df)
  
  colors <- c(
    "C>A" = "skyblue",
    "C>G" = "green",
    "T>A" = "orange",
    "T>C" = "purple",
    "T>G" = "brown",
    "C>T at CpG" = "red",
    "C>T other" = "pink"
  )
  
  ggplot(plot_df, aes(x = Mutation, y = Percentage, fill = Subtype)) +
    geom_bar(stat = "identity") +
    labs(
      title = paste0(sum(rowSums(type_occurrences_df)), " mutations"),
      x = "Mutation Type",
      y = "Percentage (%)",
      fill = "Subtype"
    ) +
    theme_classic() +
    geom_text(
      aes(label = paste0(round(Percentage, 1), "%")),
      position = position_stack(vjust = 0.5),
      size = 3,
      color = "black"
    ) +
    scale_fill_manual(values = colors) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
}

###

setwd(data_directory)
all_dirs <- list.files()
ref_genome <- "hg38"

sample_names <- all_dirs
pf_dn_grl <- GRangesList(
  setNames(
    replicate(length(sample_names), GRanges(), simplify = FALSE),
    sample_names
))

pm_dn_grl <- GRangesList(
  setNames(
    replicate(length(sample_names), GRanges(), simplify = FALSE),
    sample_names
))

for (i in 1:length(all_dirs)){
  setwd(all_dirs[i])
  all_data <- list.files()
  load(file = all_data[grep('de_novo_variant_granges', all_data)])
  pf_dn_grl[[i]] <- de_novos_pf
  pm_dn_grl[[i]] <- de_novos_pm
  
  setwd(data_directory)
}

genome(pf_dn_grl) <- ref_genome
chromosomes <- paste0('chr', c(1:22,'X', 'Y')) 
seqlevels(pf_dn_grl, pruning.mode = 'tidy') <- chromosomes

pf_dn_grl <- endoapply(pf_dn_grl, add_ref_alt)
pf_dn_grl_snvs <- get_mut_type(pf_dn_grl, type = "snv")
type_occurrences_pf <- mut_type_occurrences(pf_dn_grl_snvs, ref_genome)
type_occurrences_pf$`C>T` <- NULL
p1 <- plot_mutation_aggregate_percentages(type_occurrences_pf)


genome(pm_dn_grl) <- ref_genome
chromosomes <- paste0('chr', c(1:22,'X', 'Y')) 
seqlevels(pm_dn_grl, pruning.mode = 'tidy') <- chromosomes

pm_dn_grl <- endoapply(pm_dn_grl, add_ref_alt)
pm_dn_grl_snvs <- get_mut_type(pm_dn_grl, type = "snv")
type_occurrences_pm <- mut_type_occurrences(pm_dn_grl_snvs, ref_genome)
type_occurrences_pm$`C>T` <- NULL
p2 <- plot_mutation_aggregate_percentages(type_occurrences_pm)

