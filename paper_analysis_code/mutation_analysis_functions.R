library(dplyr)
library(tidyr)
library(ggplot2)

add_ref_alt <- function(gr) {
  # Extract the REF and ALT from the names (assuming names are in the format "chr1_3892395_C_A")
  ref_alt <- sub(".*_[^_]+_([^_]+)_([^_]+)$", "\\1,\\2", names(gr))
  ref_alt_split <- do.call(rbind, strsplit(ref_alt, ","))
  
  # Add metadata columns
  mcols(gr)$REF <- ref_alt_split[,1]
  mcols(gr)$ALT <- ref_alt_split[,2]
  
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
  
  # Handle C>T subtypes
  ct_cpg <- mutation_df$Count[mutation_df$Mutation == "C>T at CpG"]
  ct_other <- mutation_df$Count[mutation_df$Mutation == "C>T other"]
  
  # Remove C>T subtypes from the main data frame
  mutation_df <- mutation_df %>%
    filter(!Mutation %in% c("C>T at CpG", "C>T other"))
  
  # Create a new data frame for C>T with its subtypes
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
  
  # Create the bar plot
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

