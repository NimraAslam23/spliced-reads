setwd("/Users/nimraaslam/Documents/GitHub/spliced-reads/lcm")
getwd()
library(tidyverse)
library(ggplot2)
library(tidyr)
library(dplyr)

# Data import -------------------------------------------------------------

spliced_reads_orig <- read.csv("spliced_read_from_lcm.csv")
spliced_reads <- spliced_reads_orig

# Looking at data ---------------------------------------------------------

spliced_reads |> 
  distinct(sample_name) #21 unique samples

spliced_reads |> 
  distinct(junction_name) #770 splice junctions

group_column <- spliced_reads |> 
  select(junction_name, sample_name, n_spliced_reads) |> 
  complete(junction_name, sample_name,fill=list(n_spliced_reads = 0)) |> 
  mutate(disease = ifelse(grepl("ALS", sample_name),
                          "ALS",
                          "control"))

# Density plot - number of spliced reads in ALS and controls --------------

group_column |> 
  ggplot(aes(x = n_spliced_reads,
             fill = disease)) +
  stat_density() + 
  scale_x_continuous(trans = scales::pseudo_log_trans()) +
  labs(
    title = "Number of spliced reads in ALS and controls",
    x = "Number of spliced reads",
    color = "Group",
  )

# Table of the average count of each junction in each disease -------------

mean_per_junction <- group_column |> 
  group_by(junction_name,disease) |> 
  mutate(mean_n_spliced_reads = mean(n_spliced_reads)) |> 
  ungroup() |> 
  select(junction_name,disease,mean_n_spliced_reads) |> 
  unique() |> #mean_n_spliced_reads for each junction_name separately for control and ALS
  pivot_wider(names_from = 'disease',
              values_from = 'mean_n_spliced_reads')

# Nest data by junction_name ----------------------------------------------

group_column_nested <- group_column |> 
  group_by(junction_name) |> 
  nest()

# Run wilcoxon test on all the nested data --------------------------------

wilcox_tested_nested = group_column_nested |> 
  mutate(wc_res = map(data,~{broom::tidy(wilcox.test(.x$n_spliced_reads ~ .x$disease, exact=FALSE))})) |> 
  unnest(wc_res) |> 
  arrange(p.value)

# Find the junctions which were significant (p < 0.05) --------------------

significant_junctions <- wilcox_tested_nested |> 
  filter(p.value < 0.05) # 33 significant junctions

# Fraction of significant junctions expressed more in ALS and cont --------

mean_per_junction_sig <- mean_per_junction |> 
  semi_join(significant_junctions, by=("junction_name")) |> 
  mutate(higher_value = ifelse(ALS>control, "ALS", "control")) |> 
  janitor::tabyl(higher_value)

# Number of significant spliced reads in ALS and controls -----------------

sig_reads_groups <- group_column |> 
  semi_join(significant_junctions, by=("junction_name")) 

ggplot(sig_reads_groups, aes(x = n_spliced_reads, color = disease, fill = disease)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(trans = scales::pseudo_log_trans()) +
  labs(
    title = "Number of significant spliced reads in ALS and controls",
    x = "Number of spliced reads",
  )

# Junctions higher in the controls vs ALS ---------------------------------

mean_per_junction |> 
  semi_join(significant_junctions, by=("junction_name")) |> 
  mutate(higher_value = ifelse(ALS>control,"ALS","control")) |> 
  separate(junction_name, sep = '\\|',into = c("gene","junc_cat","n_datasets_junction_found"),convert = TRUE) |>
  janitor::tabyl(junc_cat,higher_value)

mean_per_junction |> 
  semi_join(significant_junctions, by=("junction_name")) |> 
  mutate(higher_value = ifelse(ALS>control,"ALS","control")) |> 
  separate(junction_name, sep = '\\|',into = c("gene","junc_cat","n_datasets_junction_found"),convert = TRUE) |> 
  ggplot(aes(x = junc_cat, fill = higher_value)) +
  geom_bar() +
  labs(
    title = "Number of junctions significantly higher in ALS and controls",
    x = "Junction category",
    y = "Number of junctions",
    color = "Group",
  )

# Mean 'n_datasets_junction_found' by which tissue it was higher in -------

mean_per_junction |> 
  semi_join(significant_junctions, by=("junction_name")) |> 
  mutate(higher_value = ifelse(ALS>control,"ALS","control")) |> group_by(higher_value) |> 
  separate(junction_name, sep = '\\|',into = c("gene","junc_cat","n_datasets_junction_found"),convert = TRUE) |> 
  summarize(mean_n_datasets = mean(n_datasets_junction_found))
