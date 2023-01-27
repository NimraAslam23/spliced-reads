getwd()
library(tidyverse)
library(ggplot2)
library(tidyr)
library(dplyr)

# Data import -------------------------------------------------------------

spliced_reads_orig <- read.csv("spliced_read_from_lcm.csv")
View(spliced_reads_orig)
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

View(group_column)



group_column |>  filter(grepl("STMN2",junction_name)) |> 
  ggplot(aes(x = n_spliced_reads,
             fill = disease)) +
  stat_density() + 
  scale_x_continuous(trans = scales::pseudo_log_trans())   

STMN2 <- group_column |>  filter(grepl("STMN2",junction_name)) 

k <- wilcox.test(n_spliced_reads ~ disease, STMN2) |> 
  broom::tidy()

# Nest data by junction_name

group_column_nested <- group_column |> 
  group_by(junction_name) |> 
  nest()

View(group_column_nested)

# Run wilcoxon test on all the nested data

mean_per_junction <- group_column_nested |> 
  mutate(mean_n_spliced_reads = map_dbl(data, ~{mean(.x$n_spliced_reads)}))

View(mean_per_junction)

wilcox.test(mean_n_spliced_reads ~ ~{(.x$disease)}, mean_per_junction)
### How to get the 'disease' part?












