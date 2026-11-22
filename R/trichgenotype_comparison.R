#****************************************************
#*This script was written by Purushoth Puvanendran and his group members and edit suggestions were proposed by Hellen Yichen Shen.
#*Question: 
#*If we compare the retrotransposon and DNA transposons content between a T. spiralis isolate and a T. nelsoni isolate, differences may be elucidated key classes that contribute to pathogencity mechanisms of interest. 
#*
#*Therefore, how do transposon content vary between T. spiralis isolate and a T. nelsoni isolate?
#*
#****************************************************

## _ Packages used -------

library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())

# Reading data and creating data frames 

setwd("./R/")

df_spiralis <- read_csv(file = "../data/trichinella_spiralis_te.csv") 
df_nelsoni <- read_csv(file = "../data/trichinella_nelsoni_te.csv") 

# Identify all the unique matching classes and visually inspect they are the same for both species 

sort(unique(df_nelsoni$matching_class))
sort(unique(df_spiralis$matching_class))

# Create vectors with DNA transposon and retrotransposon matching classes as identified by Copilot

ls_dna_transposon = c("DNA/hAT","DNA/hAT-Ac","DNA/MULE-MuDR","DNA/PIF-Harbinger","DNA/TcMar-ISRm11","DNA/TcMar-Mariner","DNA/TcMar-Pogo","DNA/TcMar-Tc1","DNA/TcMar-Tc4")

ls_retrotransposon = c("LINE/CR1","LINE/RTE-BovB","LTR/Gypsy","LTR/Ngaro","LTR/Pao","PLE","PLE/Naiad","SINE/tRNA-Deu-L2","SINE/tRNA-RTE")


## Edit 1: Reduce code redundancy using a function
# Create a function to summarize TE content for each species

summarise_te <- function(df, species_name) {
  
  df |>
    mutate(
      te_type = case_when(
        matching_class %in% ls_dna_transposon ~ "DNA transposon",
        matching_class %in% ls_retrotransposon ~ "Retrotransposon",
        TRUE ~ "Other"
      )
    ) |>
    filter(te_type != "Other") |>
    group_by(te_type) |>
    summarise(
      n_copies = n(),
      total_bp = sum(qry_width),
      .groups = "drop"
    ) |>
    mutate(
      pct_copies = 100 * n_copies / sum(n_copies),
      pct_bp = 100 * total_bp / sum(total_bp),
      species = species_name
    )
}

# Apply the function to both species

df_spiralis_summary <- summarise_te(df_spiralis, "spiralis")
df_nelsoni_summary <- summarise_te(df_nelsoni, "nelsoni")

# Combine the summary tables

df_te_summary <- bind_rows(
  df_spiralis_summary,
  df_nelsoni_summary
)




## Edit 2: Compare TE copy abundance and sequence length
# Compare TE copy abundance and sequence length

df_te_plot <- df_te_summary |>
  select(species, te_type, pct_copies, pct_bp) |>
  pivot_longer(
    cols = c(pct_copies, pct_bp),
    names_to = "measurement",
    values_to = "percentage"
  )

# Plot both measurements side by side

ggplot(df_te_plot,
       aes(x = species, y = percentage, fill = te_type)) +
  geom_col() +
  facet_wrap(~measurement) +
  labs(
    x = "Species",
    y = "Percentage (%)",
    fill = "TE type",
    title = "TE Copy Abundance and Sequence Length"
  ) +
  theme_minimal()



# Perform a Chi-squared statistical test 

counts <- matrix(df_te_summary$n_copies, nrow = 2, byrow = TRUE)
chisq.test(counts)

# Results: X-squared = 0.24545, df = 1, p-value = 0.6203
# Conclusion: T. spiralis and T. nelsoni have statistically indistinguishable proportions of DNA transposons and retrotransposons.
#******************************************************************************
## Figure 3: Retrotransposon sequence divergence in T. spiralis vs. T. nelsoni ---- 

# Creating a subset of retrotransposon matching classes that show visually observable differences in abundance between T. spiralis and T. nelsoni 

ls_retro_subset = c("PLE","PLE/Naiad","LTR/Pao", "LINE/CR1", "LTR/Gypsy")

# Filter the datasets to have only the retrotransposons of interest

df_spiralis_div <- df_spiralis |>
  filter(matching_class %in% ls_retro_subset) |>
  mutate(
    species = "T. spiralis"
  ) 

# Repeat for T. nelsoni 

df_nelsoni_div <- df_nelsoni |>
  filter(matching_class %in% ls_retro_subset) |>
  mutate(
    species = "T. nelsoni"
  ) 

# Combining dataframes 

df_trich_div <- rbind(df_spiralis_div, df_nelsoni_div)

count(df_trich_div, species, matching_class)

# Plotting a histogram to show repeat landscape 

ggplot(df_trich_div, aes(x = perc_div, fill = species)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.5) +
  facet_wrap(~matching_class) +
  labs(x = "% divergence from consensus", y = "No. of copies", fill = NULL,
       title = "Repeat landscape for select retrotransposons in Trichinella species") +
  theme_minimal()

# Performing a Kolmogorov-Smirnov (KS) test to see if both distributions come from the same underlying population 

## Edit 3: Automate KS tests using a loop
# Run KS tests for all selected retrotransposon classes

ls_ks_results <- list()

for (te in ls_retro_subset) {
  
  spiralis_values <- df_spiralis_div |>
    filter(matching_class == te) |>
    pull(perc_div)
  
  nelsoni_values <- df_nelsoni_div |>
    filter(matching_class == te) |>
    pull(perc_div)
  
  ks_result <- ks.test(
    spiralis_values,
    nelsoni_values
  )
  
  # Store the p-value for each TE class
  ls_ks_results[[te]] <- ks_result$p.value
}

# Create a summary table and adjust p-values

df_ks_results <- tibble(
  te_class = names(ls_ks_results),
  p_value = unlist(ls_ks_results)
) |>
  mutate(
    adjusted_p = p.adjust(
      p_value,
      method = "bonferroni"
    )
  )

print(df_ks_results)


