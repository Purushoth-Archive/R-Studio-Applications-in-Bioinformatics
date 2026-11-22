#****************************************************
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

df_spiralis <- read_csv(file = "../data/trichinella_spiralis_te.csv") 
df_nelsoni <- read_csv(file = "../data/trichinella_nelsoni_te.csv") 

# Identify all the unique matching classes and visually inspect they are the same for both species 

sort(unique(df_nelsoni$matching_class))
sort(unique(df_spiralis$matching_class))

# Create vectors with DNA transposon and retrotransposon matching classes as identified by Copilot

ls_dna_transposon = c("DNA/hAT","DNA/hAT-Ac","DNA/MULE-MuDR","DNA/PIF-Harbinger","DNA/TcMar-ISRm11","DNA/TcMar-Mariner","DNA/TcMar-Pogo","DNA/TcMar-Tc1","DNA/TcMar-Tc4")

ls_retrotransposon = c("LINE/CR1","LINE/RTE-BovB","LTR/Gypsy","LTR/Ngaro","LTR/Pao","PLE","PLE/Naiad","SINE/tRNA-Deu-L2","SINE/tRNA-RTE")

# Create a new column to identify which query is a DNA transposon vs. retrotransposon and create a summary table with the counts for each type of TE, % of what share of all the transposon insertions belongs to each TE type, and % of what share of the transposon-derived sequence belongs to each TE type

df_spiralis_summary <- df_spiralis |>
  mutate(
    te_type = case_when(
      matching_class %in% ls_dna_transposon ~ "DNA transposon",
      matching_class %in% ls_retrotransposon ~ "Retrotransposon",
      TRUE ~ "Other"
    )
  ) |>
  filter(te_type != "Other") |>
  group_by(te_type) |>
  summarise(n_copies = n(), total_bp = sum(qry_width)) |>
  mutate(
    pct_copies = 100 * n_copies / sum(n_copies),
    pct_bp = 100 * total_bp / sum(total_bp),
    species = "spiralis"
  )

# Repeating above coding steps for T. nelsoni

df_nelsoni_summary <- df_nelsoni |>
  mutate(
    te_type = case_when(
      matching_class %in% ls_dna_transposon ~ "DNA transposon",
      matching_class %in% ls_retrotransposon ~ "Retrotransposon",
      TRUE ~ "Other"
    )
  ) |>
  filter(te_type != "Other") |>
  group_by(te_type) |>
  summarise(n_copies = n(), total_bp = sum(qry_width)) |>
  mutate(
    pct_copies = 100 * n_copies / sum(n_copies),
    pct_bp = 100 * total_bp / sum(total_bp),
    species = "nelsoni"
  )

# Create a new data frame which combine the data in the summary data frames created above

df_te_summary <- rbind(df_spiralis_summary, df_nelsoni_summary)

# Plotting DNA transposons vs. retrotransposons in T. spiralis and T. nelsoni

ggplot(df_te_summary, aes(x = species, y = pct_copies, fill = te_type)) +
  geom_col() +
  labs(x = NULL, y = "% of TE copies", fill = NULL,
       title = "DNA transposons vs. Retrotransposons in Trichinella species") +
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

# LINE/CR1 class 

spiralis_cr1 <- df_spiralis_div |>
  filter(matching_class == "LINE/CR1")

nelsoni_cr1 <- df_nelsoni_div |>
  filter(matching_class == "LINE/CR1")

ks_cr1 <- ks.test(spiralis_cr1$perc_div, nelsoni_cr1$perc_div)
ks_cr1
rm(spiralis_cr1, nelsoni_cr1)

# LTR/Pao class

spiralis_pao <- df_spiralis_div |>
  filter(matching_class == "LTR/Pao")

nelsoni_pao <- df_nelsoni_div |>
  filter(matching_class == "LTR/Pao")

ks_pao <- ks.test(spiralis_pao$perc_div, nelsoni_pao$perc_div)
ks_pao
rm(spiralis_pao, nelsoni_pao)

# PLE class

spiralis_ple <- df_spiralis_div |>
  filter(matching_class == "PLE")

nelsoni_ple <- df_nelsoni_div |>
  filter(matching_class == "PLE")

ks_ple <- ks.test(spiralis_ple$perc_div, nelsoni_ple$perc_div)
ks_ple
rm(spiralis_ple, nelsoni_ple)

# PLE/Naiad class

spiralis_nai <- df_spiralis_div |>
  filter(matching_class == "PLE/Naiad")

nelsoni_nai <- df_nelsoni_div |>
  filter(matching_class == "PLE/Naiad")

ks_nai <- ks.test(spiralis_nai$perc_div, nelsoni_nai$perc_div)
rm(spiralis_nai, nelsoni_nai)

# LTR/Gypsy class

spiralis_gypsy <- df_spiralis_div |>
  filter(matching_class == "LTR/Gypsy")

nelsoni_gypsy <- df_nelsoni_div |>
  filter(matching_class == "LTR/Gypsy")

ks_gypsy <- ks.test(spiralis_gypsy$perc_div, nelsoni_gypsy$perc_div)
ks_gypsy
rm(spiralis_gypsy, nelsoni_gypsy)

# Adjusting p-values to ensure there are no false positives by multiplying by the number of KS tests performed (4)

ls_ks_pvalues <- c("LINE/CR1" = ks_cr1$p.value, "LTR/Pao" = ks_pao$p.value, "PLE" = ks_ple$p.value, "PLE/Naiad" = ks_nai$p.value, "LTR/Gypsy" = ks_gypsy$p.value)

p.adjust(ls_ks_pvalues, method = "bonferroni")

rm(ks_cr1, ks_nai, ks_pao, ks_ple)
